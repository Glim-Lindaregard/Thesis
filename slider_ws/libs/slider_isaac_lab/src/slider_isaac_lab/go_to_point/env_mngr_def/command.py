from isaaclab.utils import configclass
import isaaclab.envs.mdp as mdp


@configclass
class CommandsCfg:
    """Command terms for the MDP."""
    null = mdp.NullCommandCfg()