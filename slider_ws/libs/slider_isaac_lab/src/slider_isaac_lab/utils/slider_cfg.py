from dataclasses import dataclass
import torch


@dataclass
class SliderCfg:
    """
    Data class for slider robot parameters
    """
    force_max: torch.Tensor = torch.tensor(0.7, dtype=torch.float) # N
    torque_max: torch.Tensor = torch.tensor(0.0098, dtype=torch.float) # Nm
    mass: torch.Tensor = torch.tensor(4.528, dtype=torch.float) # Kg
