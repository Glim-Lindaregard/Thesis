from isaaclab.utils import configclass
import isaaclab.envs.mdp as mdp
from isaaclab.managers import SceneEntityCfg
from isaaclab.managers import TerminationTermCfg as DoneTerm


@configclass
class TerminationsCfg:
    """Termination terms for the MDP."""

    time_out = DoneTerm(func=mdp.time_out, time_out=True)

    slider_out_of_bounds = DoneTerm(
        func=mdp.illegal_contact,
        params={"threshold": 0.001, "sensor_cfg": SceneEntityCfg("contact_forces")},
    )