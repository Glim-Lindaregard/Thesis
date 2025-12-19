from .env_mngr_def.action import ActionsCfg
from .env_mngr_def.scene import SliderSceneCfg
from .env_mngr_def.observation import ObservationsCfg
from .env_mngr_def.events import EventCfg
from .env_mngr_def.termination import TerminationsCfg
from .env_mngr_def.command import CommandsCfg
from .env_mngr_def.curriculum import CurriculumCfg
from .env_mngr_def.reward import RewardsCfg
from .env_mngr_def.slider_go_to_point_env import SliderEnvGoToPoint, SliderEnvCfg
# from .env_dir_def.slider_go_to_point_env import SliderGoToPointEnv
import gymnasium as gym

gym.register(
    id="Slider-GoToPoint-v0",
    entry_point="isaaclab.envs:ManagerBasedRLEnv",
    disable_env_checker=True,
    kwargs={
        "env_cfg_entry_point": "slider_isaac_lab.go_to_point.env_mngr_def.slider_go_to_point_env:SliderEnvCfg",
        # "env_cfg_entry_point": "slider_isaac_lab.go_to_point.env_dir_def.slider_go_to_point_env:SliderGoToPointEnvCfg",
        # "env_cfg_entry_point": f"{__name__}.cartpole_env_cfg:CartpoleEnvCfg",
        "sb3_cfg_entry_point": "slider_isaac_lab.go_to_point.configs:sb3_ppo_cfg.yaml",
        # Optional training presets:
        # "skrl_cfg_entry_point": "slider_isaac_lab.agents:skrl_ppo_cfg.yaml",
        # "rsl_rl_cfg_entry_point": "slider_isaac_lab.agents.rsl_rl_ppo_cfg:PPORunnerCfg",
        # "rl_games_cfg_entry_point": "slider_isaac_lab.agents:rl_games_ppo_cfg.yaml",
    },
)
