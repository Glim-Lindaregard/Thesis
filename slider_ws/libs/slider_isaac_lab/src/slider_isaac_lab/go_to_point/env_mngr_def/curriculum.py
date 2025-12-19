from isaaclab.utils import configclass
from isaaclab.managers import CurriculumTermCfg as CurrTerm
import isaaclab.envs.mdp as mdp
import torch
from torch.distributions import Beta

def set_new_omega(env, env_ids, old_value):
    total_steps = env.common_step_counter * env.num_envs
    
    if total_steps <= 10_000_000:
        return 0.0
    else:
        new_omega = 0.003 * total_steps / 1_000_000

        alpha = 3.0
        beta = 0.9 
        r = Beta(alpha, beta).sample().item()
        new_omega = r * new_omega
        new_omega = min(new_omega, 0.2)
        return new_omega

@configclass
class CurriculumCfg:
    """Configuration for the curriculum."""
    
    new_omega_lissajou = CurrTerm(
        func=mdp.modify_env_param,
        params={
            "address": "event_manager.cfg.set_target_velocity_per_step.params.omega",
            "modify_fn": set_new_omega,
            "modify_params": {
            },
        },
    )
