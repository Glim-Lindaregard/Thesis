from isaaclab.utils import configclass
from isaaclab.managers import EventTermCfg as EventTerm
import isaaclab.envs.mdp as mdp
from isaaclab.managers import SceneEntityCfg
import math, torch
from isaaclab.envs import ManagerBasedEnv
from typing import Tuple

def sample_lissajous_params(
    env: ManagerBasedEnv,
    env_ids: torch.Tensor,
    A_range: Tuple,
    B_range: Tuple,
    a_range: Tuple,
    b_range: Tuple,
    delta_range: Tuple,
):
    """
    Sample Lissajous params once on reset
    and store in buffer for use during episode execution for each environment 
    """    
    env.lissajou_param_buffer = {}
    # physical time elapsed in environment
    env.lissajou_param_buffer["lissa_t"] = torch.zeros(env.num_envs, device=env.device)
    
    # Lissajou trajectory defined as: 
    # x = A*sin(a*t+delta)
    # y = B*sin(b*t)
    env.lissajou_param_buffer["lissa_A"] = torch.empty(env.num_envs, device=env.device).uniform_(A_range[0], A_range[1])
    env.lissajou_param_buffer["lissa_B"] = torch.empty(env.num_envs, device=env.device).uniform_(B_range[0], B_range[1])
    env.lissajou_param_buffer["lissa_a"] = torch.randint(int(a_range[0]), int(a_range[1]) + 1, (env.num_envs,), device=env.device)
    env.lissajou_param_buffer["lissa_b"] = torch.randint(int(b_range[0]), int(b_range[1]) + 1, (env.num_envs,), device=env.device)
    env.lissajou_param_buffer["lissa_delta"] = torch.empty(env.num_envs, device=env.device).uniform_(delta_range[0], delta_range[1])


def place_on_lissajous_start(
    env: ManagerBasedEnv,
    env_ids : torch.Tensor,
    asset_cfg: SceneEntityCfg,
):
    """
    One-time pose placement on the Lissajous curve at t=0
    """
    obj = env.scene[asset_cfg.name]

    A_t = env.lissajou_param_buffer["lissa_A"]
    B_t = env.lissajou_param_buffer["lissa_B"]
    d_t = env.lissajou_param_buffer["lissa_delta"]

    x0 = A_t * torch.sin(d_t)
    y0 = torch.zeros_like(B_t)

    x0 = x0[env_ids]
    y0 = y0[env_ids]
    z0 = torch.full_like(x0, 0.1)

    pos = torch.stack([x0, y0, z0], dim=-1) + env.scene.env_origins[env_ids]
    quat = torch.tensor([1.0, 0.0, 0.0, 0.0], device=env.device).expand(len(env_ids), 4)

    pose = torch.cat([pos, quat], dim=1)
    obj.write_root_pose_to_sim(pose, env_ids=env_ids)
    env.lissajou_param_buffer["lissa_t"][env_ids] = 0.0




def set_lissajous_velocity(
    env: ManagerBasedEnv,
    env_ids: torch.Tensor,
    asset_cfg: SceneEntityCfg,
    omega: float,
):
    obj = env.scene[asset_cfg.name]
    
    env.lissajou_param_buffer["lissa_t"][env_ids] += env.physics_dt
    t_all = env.lissajou_param_buffer["lissa_t"]

    theta_all = omega * t_all
    
    A_t = env.lissajou_param_buffer["lissa_A"]
    B_t = env.lissajou_param_buffer["lissa_B"]
    a_t = env.lissajou_param_buffer["lissa_a"]
    b_t = env.lissajou_param_buffer["lissa_b"]
    d_t = env.lissajou_param_buffer["lissa_delta"]

    A_t = A_t[env_ids]
    B_t = B_t[env_ids]
    a_t = a_t[env_ids]
    b_t = b_t[env_ids]
    d_t = d_t[env_ids]
    theta = theta_all[env_ids]
    m = len(env_ids)
    

    vx = A_t * a_t * omega * torch.cos(a_t * theta + d_t)
    vy = B_t * b_t * omega * torch.cos(b_t * theta)
    
    vz = torch.zeros_like(vx)
    lin = torch.stack([vx, vy, vz], dim=-1)
    ang = torch.zeros((m, 3), device=env.device)
    vel = torch.cat([lin, ang], dim=-1)
    
    obj.write_root_velocity_to_sim(vel, env_ids=env_ids)



@configclass
class EventCfg:
    """Configuration for events."""
        
    reset_slider_to_random = EventTerm(
        func=mdp.reset_root_state_uniform,
        mode="reset",
        params={
            "asset_cfg": SceneEntityCfg("slider"),
            "pose_range": {
                "x": (-1.7, 1.7), "y": (-1.7, 1.7), "z": (0.1, 0.1),
                "roll":(0.0, 0.0), "pitch": (0.0, 0.0), "yaw":(0.0, 0.0)
            },
            "velocity_range" : {
                "x": (0, 0), "y": (0, 0), "z": (0.0, 0.0),
                "roll":(0.0, 0.0), "pitch": (0.0, 0.0), "yaw":(0.0, 0.0)
            }
        }
    )

    set_target_trajectory_params_at_reset = EventTerm(
        func=sample_lissajous_params,
        mode="reset",
        params={
            "A_range": (0.5, 1.9),
            "B_range": (0.5, 1.9),
            "a_range": (1.0, 3.0),
            "b_range": (1.0, 3.0),
            "delta_range": (0.0, 2.0 * math.pi),
        },
    )

    set_target_pose_at_reset = EventTerm(
        func=place_on_lissajous_start,
        mode="reset",
        params={
            "asset_cfg": SceneEntityCfg("target_sphere"),
        },
    )

    set_target_velocity_per_step = EventTerm(
        func=set_lissajous_velocity,
        mode="interval",
        interval_range_s=(0.0, 0.0),
        params={
            "asset_cfg": SceneEntityCfg("target_sphere"),
            "omega": 0.0,  # adjust via curriculum learning
        },
    )

