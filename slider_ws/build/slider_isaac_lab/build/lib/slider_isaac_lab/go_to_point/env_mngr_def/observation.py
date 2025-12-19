from isaaclab.envs import ManagerBasedEnv
from isaaclab.managers import SceneEntityCfg
from isaaclab.assets import RigidObject
from isaaclab.managers import ObservationGroupCfg as ObsGroup
from isaaclab.managers import ObservationTermCfg as ObsTerm
from isaaclab.utils import configclass

from slider_isaac_lab.utils.get_admissible_control import get_admissible_action_bounds
from slider_isaac_lab.utils.slider_cfg import SliderCfg
import torch


from isaaclab.utils import math as il_math  


def quat_to_rpy_wxyz(q):
    w, x, y, z = q.unbind(-1)
    roll = torch.atan2(2*(w*x + y*z), 1 - 2*(x*x + y*y))
    pitch = torch.asin(2*(w*y - z*x))
    yaw = torch.atan2(2*(w*z + x*y), 1 - 2*(y*y + z*z))
    return torch.stack((roll, pitch, yaw), dim=-1)


def get_admissible_action_bounds_wrapper(env: ManagerBasedEnv, 
                    slider_cfg: SceneEntityCfg,
                    target_cfg: SceneEntityCfg,
                    ) -> torch.Tensor:  #(B,2,2)
    """
    Get admissible action bounds based on analytical formulas 
    for fx, fy world frame forces
    """
    
    slider: RigidObject = env.scene[slider_cfg.name]
    force_max = SliderCfg.force_max

    base = torch.tensor([[-force_max, force_max],[-force_max, force_max]],
                        dtype=torch.float,
                        device=slider.data.root_lin_vel_w.device)   

    limit = base.unsqueeze(0).expand(env.scene.env_origins[:, :2].shape[0], -1, -1)
    
    pos = slider.data.root_pos_w[:, :2] - env.scene.env_origins[:, :2]       
    vel = slider.data.root_lin_vel_w[:, :2]  

    B = pos.shape[0]
    x = torch.zeros(B, 6, 
                    device=pos.device, 
                    dtype=pos.dtype)  
    
    q = slider.data.root_state_w[:, 3:7]   
    rpy = quat_to_rpy_wxyz(q)

    x[:, 0:2] = pos     
    x[:, 3:5] = vel     
    x[:, 2] = rpy[:, 2]

    box_position_xy_bounds = torch.as_tensor(
        [[-2, 2],
        [-2, 2]],
        device=x.device,
        dtype=x.dtype,
    )
    
    limit = get_admissible_action_bounds(
        x=x, 
        x_bounds=box_position_xy_bounds, 
        )
    
    return limit


def get_observations(
    env: ManagerBasedEnv,
    slider_cfg: SceneEntityCfg,
    target_cfg: SceneEntityCfg,
    ) -> torch.Tensor: #(B,8)
    
    """
    Get Policy observations:
    - True position and velocity of robot
    - Reference position velocity relative to target point
    
    """

    slider: RigidObject = env.scene[slider_cfg.name]
    target: RigidObject = env.scene[target_cfg.name]

    pos = slider.data.root_pos_w[:, :2] \
        - env.scene.env_origins[:, :2]     

    vel = slider.data.root_lin_vel_w[:, :2]

    target_pos = target.data.root_pos_w[:, :2] \
        - env.scene.env_origins[:, :2]

    err_pos = pos - target_pos 
    err_vel = vel - target.data.root_lin_vel_w[:, :2]

    return torch.cat([pos, vel, err_pos, err_vel], dim=-1)


@configclass
class ObservationsCfg:
    """Observation specifications for the MDP."""

    @configclass
    class PolicyCfg(ObsGroup):
        """Observations for policy group."""        
        
        observations = ObsTerm(
            func=get_observations,
            params={
                "slider_cfg": SceneEntityCfg("slider"),
                "target_cfg": SceneEntityCfg("target_sphere"),
            },
        )
    
        action_bounds = ObsTerm(
            func=get_admissible_action_bounds_wrapper,
            params={
                "slider_cfg": SceneEntityCfg("slider"),
                "target_cfg": SceneEntityCfg("target_sphere"),
            },
        )
        def __post_init__(self):
            self.enable_corruption = False
            self.concatenate_terms = False

    policy: PolicyCfg = PolicyCfg()