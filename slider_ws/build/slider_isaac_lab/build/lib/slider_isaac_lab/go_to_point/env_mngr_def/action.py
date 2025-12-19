from __future__ import annotations
from isaaclab.managers import ActionTerm, ActionTermCfg
from isaaclab.assets import RigidObject
from isaaclab.envs import ManagerBasedEnv
from isaaclab.utils import configclass
import torch 



class SliderActionTerm(ActionTerm):
    _asset: RigidObject
    """The rigid body asset on which the action term is applied."""
    
    def __init__(self, cfg: SliderActionTermCfg, env: ManagerBasedEnv):
        super().__init__(cfg, env)
        
        self._raw_actions = torch.zeros(env.num_envs, 2, device=self.device) 
        self._processed_actions = torch.zeros(env.num_envs, 6, device=self.device)
        
    @property
    def total_action_dim(self) -> int:
        return self._raw_actions.shape[1]

    @property
    def action_dim(self) -> int:
        return self._raw_actions.shape[1]
    
    @property
    def raw_actions(self) -> torch.Tensor:
        return self._raw_actions
    
    @property
    def processed_actions(self) -> torch.Tensor:
        return self._processed_actions
    
    def process_actions(self, actions: torch.Tensor):
        self._raw_actions = actions
        self._processed_actions[:,:2] = self._raw_actions[:]
        
    def apply_actions(self):
        self._asset.set_external_force_and_torque(
            forces=self._processed_actions[:,:3].unsqueeze(1), 
            torques=self._processed_actions[:,3:].unsqueeze(1),
            is_global=True)
        self._asset.write_data_to_sim()

@configclass
class SliderActionTermCfg(ActionTermCfg):
    class_type: type = SliderActionTerm
    """The class corresponding to the action term."""


@configclass
class ActionsCfg:
    """Action specifications for the MDP"""
    joint_pos = SliderActionTermCfg(asset_name="slider")