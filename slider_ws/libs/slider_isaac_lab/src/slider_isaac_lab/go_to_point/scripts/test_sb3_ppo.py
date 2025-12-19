import argparse
from isaaclab.app import AppLauncher
import gym
import torch
from stable_baselines3 import PPO


# --- Args ---
parser = argparse.ArgumentParser()
parser.add_argument("--num_envs", type=int, default=1)
parser.add_argument("--policy_path", type=str, required=True, help="Path to PPO policy .zip")
AppLauncher.add_app_launcher_args(parser)
args_cli = parser.parse_args()
args_cli.livestream = 2

# --- Launch app ---
app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app


from isaaclab.envs import ManagerBasedRLEnvCfg, ManagerBasedRLEnv
from isaaclab.utils import configclass
from isaaclab.sim import SimulationCfg, PhysxCfg

import slider_isaac_lab.go_to_point


# --- Main ---
def main():
    env_cfg = parse_env_cfg(
        args_cli.task, device=args_cli.device, num_envs=args_cli.num_envs, use_fabric=not args_cli.disable_fabric
    )
    
    env = gym.make(args_cli.task, cfg=env_cfg)

    model = PPO.load(args_cli.policy_path)

    obs, _ = env.reset()
    count = 0
    while simulation_app.is_running():
        with torch.inference_mode():
            action, _ = model.predict(obs, deterministic=True)
            obs, rew, terminated, truncated, info = env.step(action)
            count += 1

    env.close()


if __name__ == "__main__":
    main()
    simulation_app.close()
