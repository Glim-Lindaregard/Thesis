"""
Open-Loop Test of go-to-point environemnt 


"""

import argparse
import os
from isaaclab.app import AppLauncher

parser = argparse.ArgumentParser(description="Main script for running MDP")
parser.add_argument("--num_envs", type=int, default=1, help="Number of environments to spawn.")

# append AppLauncher cli args
AppLauncher.add_app_launcher_args(parser)
args_cli = parser.parse_args()
args_cli.livestream = 2
# public_ip = os.environ.get("PUBLIC_IP", "0.0.0.0")  # set PUBLIC_IP in your shell first

app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app
# simulation_app.set_setting("/app/livestream/allowAllHosts", True)
# simulation_app.set_setting("/app/livestream/addr", "0.0.0.0")
# simulation_app.set_setting("/app/livestream/publicEndpointAddress", public_ip)
# simulation_app.set_setting("/app/livestream/port", 49100)  # int, not string


"""Rest everything follows."""

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

    scene: SliderSceneCfg = SliderSceneCfg(num_envs=args_cli.num_envs,
                                           env_spacing=5.0, 
                                           replicate_physics=False) 
    
    observations: ObservationsCfg = ObservationsCfg()
    actions: ActionsCfg = ActionsCfg()
    events: EventCfg = EventCfg()
    
    sim: SimulationCfg = SimulationCfg(physics_prim_path="/physicsScene", 
                                       dt=0.1, 
                                       render_interval=1, 
                                       gravity=(0.0, 0.0, 0.0), 
                                       device="cuda",
                                       physx=PhysxCfg())
    
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
        self.episode_length_s = 10  


def main():
    """Main function."""

    # setup base environment
    env = ManagerBasedRLEnv(cfg=SliderEnvCfg())

    action = torch.zeros(env.num_envs, 3, device=env.device)
    
    for i in range(env.num_envs):
        action[i, :]= torch.tensor([1,0,0], device=env.device)

    count = 0
    while simulation_app.is_running():
        with torch.inference_mode():
            
            # reset
            if count % 2048 == 0:
                # count = 0
                print(count)
                # obs, _ = env.reset(seed=None)
                # print("[INFO]: Resetting environment...")

            obs, rew, terminated, truncated, info = env.step(action)            
        

            count += 1*args_cli.num_envs

    env.close()


if __name__ == "__main__":
    main()
    simulation_app.close()
