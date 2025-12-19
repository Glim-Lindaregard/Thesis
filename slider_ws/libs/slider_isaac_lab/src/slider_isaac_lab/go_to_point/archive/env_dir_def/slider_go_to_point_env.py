# Copyright (c) 2022-2025, The Isaac Lab Project Developers (https://github.com/isaac-sim/IsaacLab/blob/main/CONTRIBUTORS.md).
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations

import math
import torch
from collections.abc import Sequence

from isaaclab_assets.robots.cartpole import CARTPOLE_CFG

import isaaclab.sim as sim_utils
from isaaclab.assets import Articulation, ArticulationCfg
from isaaclab.envs import DirectRLEnv, DirectRLEnvCfg
from isaaclab.scene import InteractiveSceneCfg
from isaaclab.sim import SimulationCfg
from isaaclab.sim.spawners.from_files import GroundPlaneCfg, spawn_ground_plane
from isaaclab.utils import configclass
from isaaclab.utils.math import sample_uniform

import argparse
import os
from isaaclab.app import AppLauncher

from isaaclab.utils import configclass
from isaaclab.sim import SimulationCfg, PhysxCfg


# from slider_isaac_lab.src.slider_isaac_lab.go_to_point.env_dir_def.configs import SliderEnvCfg
from slider_isaac_lab.utils.get_admissible_control import get_admissible_action_bounds

from slider_isaac_lab.utils.slider_cfg import SliderCfg
from isaaclab.assets import RigidObject
import isaaclab.envs.mdp as mdp
from slider_isaac_lab.go_to_point.env_dir_def.scene import SliderSceneCfg
from gymnasium import spaces


import torch
import matplotlib.pyplot as plt 
import numpy as np 
    
        
@configclass
class SliderGoToPointEnvCfg(DirectRLEnvCfg):
    
    # env
    decimation = 1
    episode_length_s = 20.0
    action_scale = 1  
    state_space = 0
    seed = 10
    is_finite_horizon: bool = True
    
    observation_space = {"observations":8, "action_bounds":4}
        
    action_space = 2

    # simulation
    sim: SimulationCfg = SimulationCfg(physics_prim_path="/physicsScene", 
                                       render_interval=1, 
                                       gravity=(0.0, 0.0, 0.0), 
                                       device="cuda",
                                       physx=PhysxCfg())
    # scene
    scene: SliderSceneCfg = SliderSceneCfg(num_envs=1,
                                           env_spacing=10.0, 
                                           replicate_physics=False) 
    
    
    
    

class SliderGoToPointEnv(DirectRLEnv):
    cfg: SliderGoToPointEnvCfg
    
    def __init__(self, cfg: SliderGoToPointEnvCfg, render_mode: str | None = None, **kwargs):
        super().__init__(cfg, render_mode, **kwargs)

    def _setup_scene(self):
        self.slider: RigidObject = self.scene.slider
        self.target: RigidObject = self.scene.target
        
        self.scene.clone_environments(copy_from_source=False)

    def _pre_physics_step(self, actions: torch.Tensor) -> None:
        self.actions = self.action_scale * actions.clone()

    def _apply_action(self) -> None:
        
        forces_xy = self._actions
        
        forces = torch.zeros(self.num_envs, 3, device=self.device)
        forces[:, :2] = forces_xy

        self.slider.set_external_force_and_torque(forces=forces, torques=self._zero_torques)

    def _get_observations(self) -> dict:
        
        pos = self.slider.data.root_pos_w[:, :2] - self.scene.env_origins[:, :2]     
        vel = self.slider.data.root_lin_vel_w[:, :2]
    
        target_pos = self.target.data.root_pos_w[:, :2] - self.scene.env_origins[:, :2]
        target_vel = - self.target.data.root_lin_vel_w[:, :2]
        
        err_pos = pos - target_pos 
        err_vel = vel - target_vel
        
        obs = torch.cat([pos, vel, err_pos, err_vel], dim=-1)
        
        action_bounds = self.action_bounds_f()
        
        observations = {"observations": obs, 
                        "action_bounds": action_bounds}
        return observations
        

    def _get_rewards(self) -> torch.Tensor:
        
        total_reward = compute_rewards(
            slider_pos=self.slider.data.root_pos_w[:, :2],
            target_pos=self.target.data.root_pos_w[:, :2],
            reset_terminated=self.reset_terminated)
        
        return total_reward


    def _get_dones(self):
        time_out = self.episode_length_buf >= self.max_episode_length - 1

        sensor = self.scene["contact_forces"]   
        forces = sensor.data.net_forces_w       
        max_norm = torch.linalg.vector_norm(forces, dim=-1).amax(dim=-1)
        slider_out_of_bounds = max_norm > 1e-3

        return slider_out_of_bounds, time_out


    def _reset_idx(self, env_ids: Sequence[int] | None):
        if env_ids is None:
            env_ids = torch.arange(self.num_envs, device=self.device)
        elif not isinstance(env_ids, torch.Tensor):
            env_ids = torch.as_tensor(env_ids, device=self.device, dtype=torch.long)

        # Reset buffers first
        super()._reset_idx(env_ids)
        self._actions[env_ids] = 0.0

        # Current Z & orientation preserved; sample XY around each env origin
        slider_z = self.slider.data.root_pos_w[env_ids, 2]
        target_z = self.target.data.root_pos_w[env_ids, 2]
        quat_slider = self.slider.data.root_quat_w[env_ids]
        quat_target = self.target.data.root_quat_w[env_ids]

        lo_s, hi_s = self.cfg.slider_spawn_xy_range
        lo_t, hi_t = self.cfg.target_spawn_xy_range

        # Sample positions in env frame then shift by origins
        slider_xy = torch.rand(len(env_ids), 2, device=self.device) * (hi_s - lo_s) + lo_s
        target_xy = torch.rand(len(env_ids), 2, device=self.device) * (hi_t - lo_t) + lo_t

        world_slider_xy = slider_xy + self.scene.env_origins[env_ids, :2]
        world_target_xy = target_xy + self.scene.env_origins[env_ids, :2]

        # Write root state (pos, quat, lin vel, ang vel) = 13
        slider_states = torch.zeros(len(env_ids), 13, device=self.device)
        target_states = torch.zeros(len(env_ids), 13, device=self.device)

        slider_states[:, 0:3] = torch.stack([world_slider_xy[:, 0], world_slider_xy[:, 1], slider_z], dim=-1)
        slider_states[:, 3:7] = quat_slider
        # zero velocities by default

        target_states[:, 0:3] = torch.stack([world_target_xy[:, 0], world_target_xy[:, 1], target_z], dim=-1)
        target_states[:, 3:7] = quat_target

        # Push to sim
        self.slider.write_root_state_to_sim(slider_states, env_ids)
        self.target.write_root_state_to_sim(target_states, env_ids)
        
        
        
    def action_bounds_f(self) -> torch.Tensor:  #(B,3,2)

        force_max = SliderCfg.force_max

        base = torch.tensor([[-force_max, force_max],
                            [-force_max, force_max]], 
                            dtype=torch.float,
                            device=self.slider.data.root_lin_vel_w.device)  # (2,2)

        limit = base.unsqueeze(0).expand(self.scene.env_origins[:, :2].shape[0], -1, -1)
        
        pos = self.slider.data.root_pos_w[:, :2] - self.scene.env_origins[:, :2]       
        vel = self.slider.data.root_lin_vel_w[:, :2]  

        B = pos.shape[0]
        x = torch.zeros(B, 6, 
                        device=pos.device, 
                        dtype=pos.dtype)
        x[:, 0:2] = pos
        x[:, 3:5] = vel

        x_bounds = torch.as_tensor(
            [[-2, 2],
            [-2, 2]],
            device=x.device,
            dtype=x.dtype,
        )
        
        limit = get_admissible_action_bounds(x=x,
                                    x_bounds=x_bounds)
        return limit



@torch.jit.script
def compute_rewards(
    slider_pos: torch.Tensor,
    target_pos: torch.Tensor,
    reset_terminated: torch.Tensor):
    
    rew_termination =  reset_terminated.float()

    
    diff = slider_pos - target_pos
    total_reward = 1.0 / (torch.norm(diff, dim=1) + 1e-2) + rew_termination
    
    return total_reward