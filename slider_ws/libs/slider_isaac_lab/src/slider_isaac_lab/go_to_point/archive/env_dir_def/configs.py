# Copyright (c) 2022-2025, The Isaac Lab Project Developers (https://github.com/isaac-sim/IsaacLab/blob/main/CONTRIBUTORS.md).
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations

import math
import torch
from collections.abc import Sequence

import isaaclab.sim as sim_utils
from isaaclab.assets import Articulation, ArticulationCfg
from isaaclab.envs import DirectRLEnv, DirectRLEnvCfg
from isaaclab.scene import InteractiveSceneCfg
from isaaclab.sim import SimulationCfg
from isaaclab.sim.spawners.from_files import GroundPlaneCfg, spawn_ground_plane
from isaaclab.utils import configclass
from isaaclab.utils.math import sample_uniform
from isaaclab.assets import AssetBaseCfg, RigidObjectCfg

import argparse
import os
from isaaclab.app import AppLauncher

from isaaclab.envs import ManagerBasedRLEnvCfg, ManagerBasedRLEnv
from isaaclab.utils import configclass
from isaaclab.sim import SimulationCfg, PhysxCfg

from slider_isaac_lab.go_to_point import (
    ObservationsCfg,
    SliderSceneCfg,
    ActionsCfg,
    EventCfg,
    TerminationsCfg,
    RewardsCfg,
    CommandsCfg,
)

import torch
import matplotlib.pyplot as plt 
import numpy as np 
    
    
@configclass
class SliderEnvCfg(DirectRLEnvCfg):
    
    decimation = 1
    # viewer.eye = (0.0, 0.0, 30.0)
    episode_length_s = 40  


    scene: SliderSceneCfg = SliderSceneCfg(num_envs=1,
                                           env_spacing=10.0, 
                                           replicate_physics=False) 
    
    observations: ObservationsCfg = ObservationsCfg()
    actions: ActionsCfg = ActionsCfg()
    events: EventCfg = EventCfg()
    
    sim: SimulationCfg = SimulationCfg(physics_prim_path="/physicsScene", 
                                       render_interval=1, 
                                       gravity=(0.0, 0.0, 0.0), 
                                       device="cuda",
                                       physx=PhysxCfg(),
                                       dt = 0.1)    
    
    self.slider_cfg : RigidObjectCfg = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/slider", 
        spawn=sim_utils.UsdFileCfg(usd_path=slider_usd_path, 
                                   visual_material=sim_utils.PreviewSurfaceCfg(diffuse_color=(0.1, 0.1, 0.1), 
                                                                               metallic=0.2), 
                                   mass_props=sim_utils.MassPropertiesCfg(mass=4.528),
                                   activate_contact_sensors=True, 
                                   collision_props=sim_utils.CollisionPropertiesCfg(collision_enabled=True, 
                                                                                    rest_offset=0.0,
                                                                                    contact_offset=0.1)), 
        init_state=RigidObjectCfg.InitialStateCfg(pos=(0.0, 0.0, 0.1), 
                                                  rot=(1.0, 0.0, 0.0, 0.0)),
        collision_group=0,
        debug_vis=True)
    
    
    self.target_cfg: RigidObjectCfg = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/target_sphere",
        spawn=sim_utils.SphereCfg(
            radius=0.1, 
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=False),
            collision_props=sim_utils.CollisionPropertiesCfg(collision_enabled=False, 
                                                             rest_offset=0.0,
                                                             contact_offset=0.1), 
            visual_material=sim_utils.PreviewSurfaceCfg(diffuse_color=(0.5, 0.0, 0.0), 
                                                        metallic=0.0)
            ), 
        init_state=RigidObjectCfg.InitialStateCfg(pos=(0.0, 0.0, 0.1),
                                                  rot=(1.0, 0.0, 0.0, 0.0)),
        collision_group=0,
    )
