from isaaclab.utils import configclass
import isaaclab.envs.mdp as mdp
from isaaclab.managers import RewardTermCfg as RewTerm
from isaaclab.envs import ManagerBasedEnv
from isaaclab.managers import SceneEntityCfg
from isaaclab.assets import RigidObject
from isaaclab.managers import ObservationGroupCfg as ObsGroup
from isaaclab.managers import ObservationTermCfg as ObsTerm
from isaaclab.utils import configclass
import torch



def get_reward_position_error(
    env: ManagerBasedEnv, 
    slider_cfg: SceneEntityCfg,
    target_cfg: SceneEntityCfg,
    ) -> torch.Tensor: #(1,)

    slider: RigidObject = env.scene[slider_cfg.name]
    target: RigidObject = env.scene[target_cfg.name]
    
    diff = slider.data.root_pos_w[:, :2] - target.data.root_pos_w[:, :2]

    return 1.0 / (torch.norm(diff, dim=1) + 1e-2)




@configclass
class RewardsCfg:
    """Reward terms for the MDP."""
    
    # (1) Tracking reward
    tracking = RewTerm(
        func=get_reward_position_error, 
        params={"slider_cfg": SceneEntityCfg("slider"),
                "target_cfg": SceneEntityCfg("target_sphere")},
        weight=1.0
        )
    
    # (2) Failure penalty
    terminating = RewTerm(
        func=mdp.is_terminated,
        weight=-50.0
        )
    
