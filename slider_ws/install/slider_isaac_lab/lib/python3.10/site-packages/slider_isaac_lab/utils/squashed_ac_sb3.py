"""
Custom network for PPO StableBaselines3 with tanh action squashing

- BaseFeaturesExtractor: For processing dictionary actions storing both observations and action bounds
- ActorCriticPolicy: For defining the custom AC architecture with the modified continuous action distribution 

"""

# Standard library
import math
from typing import Callable, Dict, List, Optional, Tuple, Type, Union, Any

# Third-party
import gymnasium as gym 
import numpy as np
import torch as th
from stable_baselines3 import PPO
from stable_baselines3.common.distributions import Distribution
from stable_baselines3.common.policies import ActorCriticPolicy
from stable_baselines3.common.torch_layers import BaseFeaturesExtractor
from torch import nn

LOG_STD_MIN, LOG_STD_MAX = -20, 2
EPS = 1e-6


def atanh(x: th.Tensor) -> th.Tensor:
    x = x.clamp(-1 + EPS, 1 - EPS)
    return 0.5 * (th.log1p(x) - th.log1p(-x))

class BoundedTanhGaussianDistribution(Distribution):
    """
    Base Normal N(mu, sigma) -> y = tanh(z) -> a = m + s * y
    where m = (u + l)/2, s = (u - l)/2, per sample/per action.
    """
    def __init__(self, action_dim: int):
        super().__init__()
        self.mean_actions: Optional[th.Tensor] = None  # (B, A)
        self.log_std: Optional[th.Tensor] = None       # (B, A)
        self.std: Optional[th.Tensor] = None           # (B, A)
        self.bounds: Optional[th.Tensor] = None        # (B, 2, A)
        self.action_dim = action_dim
        
    def set_bounds(self, 
                bounds: th.Tensor
                )-> None:
        """
        Returns action bounds
        """
        self.bounds = bounds.detach().contiguous()

    def proba_distribution(self, 
                        mean_actions:th.Tensor, 
                        log_std: th.Tensor
                        ) -> "BoundedTanhGaussianDistribution":
        """
        Set parameters of the distribution.

        :return: self
        """
        log_std = th.clamp(log_std, LOG_STD_MIN, LOG_STD_MAX)
        if log_std.ndim == 1:
            log_std = log_std.unsqueeze(0).expand_as(mean_actions)
        self.mean_actions = mean_actions
        self.log_std = log_std
        self.std = th.exp(log_std)
        return self
    
    
    def proba_distribution_net(self, latent_dim: int =64, log_std_init: float = 0.0
                            )->Union[nn.Module, Tuple[nn.Module, nn.Parameter]]:
        """
        Create the layers and parameters that represent the distribution.

        Subclasses must define this, but the arguments and return type vary between
        concrete classes.
        """

        mean_actions = nn.Linear(latent_dim, self.action_dim)
        log_std = nn.Parameter(th.ones(self.action_dim) * log_std_init, requires_grad=True)
        
        return mean_actions, log_std
        

    def actions_from_params(self, 
                            mean_actions: th.Tensor, 
                            log_std: th.Tensor,
                            deterministic: bool = False
                            )->Any:
        """
        Returns samples from the probability distribution
        given its parameters.

        :return: actions
        """
        self.proba_distribution(mean_actions, log_std)
        if deterministic:
            actions = self.mode()  
        else:
            actions = self.sample()
        log_prob = self.log_prob(actions)
        return actions, log_prob

    def log_prob_from_params(self, 
                            mean_actions: th.Tensor, 
                            log_std: th.Tensor
                            )-> tuple[th.Tensor, th.Tensor]:
        """
        Returns samples and the associated log probabilities
        from the probability distribution given its parameters.

        :return: actions and log prob
        """
        self.proba_distribution(mean_actions, log_std)
        a = self.sample()
        return self.log_prob(a)

    def get_mid_point_scale(self):
        """
        Returns midpoint and scale to rescale tanh-bounded
        action between upper and lower bounds
        """
        l = self.bounds[..., 0]
        u = self.bounds[..., 1]
        s = 0.5 * (u - l)
        s = th.clamp(s, min=1e-4)
        m = 0.5 * (u + l)
        return m, s

    def sample(self) -> th.Tensor:
        """
        Returns a sample from the probability distribution

        :return: the stochastic action
        """
        eps = th.randn_like(self.mean_actions)
        z = self.mean_actions + self.std * eps        
        y = th.tanh(z)
        m, s = self.get_mid_point_scale()
        a = m + s * y
        return a

    def mode(self) -> th.Tensor:
        """
        Returns the most likely action (deterministic output)
        from the probability distribution

        :return: the stochastic action
        """

        z = self.mean_actions                          
        y = th.tanh(z)
        m, s = self.get_mid_point_scale()
        return m + s * y

    def log_prob(self, 
                actions: th.Tensor
                ) -> th.Tensor:
        """
        Compute log-probability of bounded actions under a tanh-squashed Gaussian policy.

        Given:
            z ~ N(mean_actions, std^2)
            y = tanh(z)
            a = m + s * y
            m = (u + l)/2
            s = (u - l)/2

        Then by the change of variables formula:
            log p(a) = log p(z)
                    - log|dy/dz|          # Tanh squashing correction
                    - log|da/dy|          # Scaling to environment bounds
        where:
            log|dy/dz| = -2 * (z + softplus(-2z))
            log|da/dy| = log(s)
        """
        
        action_midpoint, action_scale = self.get_mid_point_scale()
        y = (actions - action_midpoint) / action_scale
        y = y.clamp(-1 + EPS, 1 - EPS)
        z = atanh(y)

        # Compute Gaussian log-probability in pre-tanh space
        inv_std = 1.0 / (self.std + EPS)
        standardized_diff = (z - self.mean_actions) * inv_std
        log_prob_pre_tanh = (
            -0.5 * standardized_diff.pow(2)
            - self.log_std
            # - 0.5 * math.log(2 * math.pi)
        ).sum(dim=-1)

        # Change-of-variable correction terms
        # softplus(x) = log(1 + exp(x)) — smooth, 
        # stable version of ReLU used for numerical stability
        log_abs_dy_dz = -2.0 * (z + th.nn.functional.softplus(-2.0 * z))
        log_abs_da_dy = th.log(action_scale)

        log_prob = (
            log_prob_pre_tanh # Gaussian prob
            - log_abs_dy_dz.sum(dim=-1) # tanh correction 
            - log_abs_da_dy.sum(dim=-1) # action scale correction 
        )

        return th.nan_to_num(log_prob, nan=-1e20, neginf=-1e20, posinf=-1e20)

    def entropy(self) -> th.Tensor:
        """
        Returns Shannon's entropy of the probability

        :return: the entropy, or None if no analytical form is known
        """
        log_std = th.clamp(self.log_std, LOG_STD_MIN, LOG_STD_MAX)
        base_ent = (0.5 + 0.5 * math.log(2 * math.pi) + log_std).sum(dim=-1)
        _, s = self.get_mid_point_scale()
        return th.nan_to_num(base_ent + th.log(s).sum(dim=-1), \
            nan=0.0, neginf=0.0, posinf=0.0)

class DictObsExtractor(BaseFeaturesExtractor):
    def __init__(self, observation_space: gym.spaces.Dict, features_dim: int = 64):
        env_space: gym.spaces.Box = observation_space["observations"]
        super().__init__(observation_space, features_dim)
        self.env_dim = int(np.prod(env_space.shape))
        self.net = nn.Sequential(
            nn.Flatten(),
            nn.Linear(self.env_dim, 128), nn.ReLU(),
            nn.Linear(128, features_dim), nn.ReLU(),
        )
        self.action_dim = int(observation_space["action_bounds"].shape[-1])
        self.last_bounds: Optional[th.Tensor] = None

    def forward(self, obs: Dict[str, th.Tensor]) -> th.Tensor:
        self.last_bounds = obs["action_bounds"].float().detach()
        return self.net(obs["observations"].float())

class CustomNetwork(nn.Module):
    def __init__(self, feature_dim:int, 
                last_layer_dim_pi:int = 64,
                last_layer_dim_vf:int = 64):
        super().__init__()
        self.latent_dim_pi = last_layer_dim_pi
        self.latent_dim_vf = last_layer_dim_vf
        self.policy_net = nn.Sequential(nn.Linear(feature_dim, 
                                                last_layer_dim_pi), 
                                        nn.ReLU())
        self.value_net  = nn.Sequential(nn.Linear(feature_dim, 
                                                last_layer_dim_vf),
                                        nn.ReLU())

    def forward(self, features: th.Tensor) -> Tuple[th.Tensor, th.Tensor]:
        return self.forward_actor(features), self.forward_critic(features)

    def forward_actor(self, features: th.Tensor) -> th.Tensor:
        return self.policy_net(features)

    def forward_critic(self, features: th.Tensor) -> th.Tensor:
        return self.value_net(features)

class CustomActorCriticPolicy(ActorCriticPolicy):
    def __init__(
        self,
        observation_space: gym.spaces.Space,
        action_space: gym.spaces.Space,
        lr_schedule: Callable[[float], float],
        net_arch: Optional[List[Union[int, Dict[str, List[int]]]]] = None,
        activation_fn: Type[nn.Module] = nn.Tanh,
        *args, **kwargs,
    ):
        super().__init__(observation_space, action_space, lr_schedule, net_arch, activation_fn, *args, **kwargs)
        self.ortho_init = False
        self.action_dist = BoundedTanhGaussianDistribution(self.action_space.shape[0])
        self.log_std = nn.Parameter(th.zeros(self.action_space.shape[0]))

    def _build_mlp_extractor(self) -> None:
        """
        Create the policy and value networks.
        Part of the layers can be shared.
        """
        self.mlp_extractor = CustomNetwork(self.features_dim)
        
    def _get_action_dist_from_latent(self, 
                                    latent_pi: th.Tensor
                                    ) -> Distribution:
        """
        Retrieve action distribution given the latent codes.

        :param latent_pi: Latent code for the actor
        :return: Action distribution
        """
        mean_actions = self.action_net(latent_pi)

        log_std = self.log_std.expand_as(mean_actions)

        bounds = getattr(self.features_extractor, "last_bounds", None)
        assert bounds is not None, "Bounds missing. Make sure observation includes {'bounds': (2, action_dim)}."

        self.action_dist.set_bounds(bounds)
        return self.action_dist.proba_distribution(mean_actions, log_std)

