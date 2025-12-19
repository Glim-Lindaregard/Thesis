"""
Manager based environment definition for slider robot:

Details:
- Objective is to reach a target point 
- Environment is a box shaped arena 
- 

"""

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
    CurriculumCfg
)

import torch
import matplotlib.pyplot as plt 
import numpy as np 

@configclass
class SliderEnvCfg(ManagerBasedRLEnvCfg):

    scene: SliderSceneCfg = SliderSceneCfg(
        num_envs=1,
        env_spacing=10.0, 
        replicate_physics=False
    ) 

    
    observations: ObservationsCfg = ObservationsCfg()
    actions: ActionsCfg = ActionsCfg()
    events: EventCfg = EventCfg()

    sim: SimulationCfg = SimulationCfg(
        physics_prim_path="/physicsScene", 
        render_interval=1, 
        gravity=(0.0, 0.0, 0.0), 
        device="cuda",
        physx=PhysxCfg()
    )
    
    # MDP settings
    curriculum: CurriculumCfg = CurriculumCfg()
    rewards: RewardsCfg = RewardsCfg()
    terminations: TerminationsCfg = TerminationsCfg()
    commands: CommandsCfg = CommandsCfg()
    

    def __post_init__(self):
        """Post initialization."""
        self.decimation = 1
        self.sim.dt = 0.1
        self.viewer.eye = (0.0, 0.0, 30.0)
        self.episode_length_s = 20  



class SliderEnvGoToPoint(ManagerBasedRLEnv):
    cfg = SliderEnvCfg()
