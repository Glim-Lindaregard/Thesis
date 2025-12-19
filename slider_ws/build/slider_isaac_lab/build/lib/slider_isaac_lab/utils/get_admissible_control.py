import torch
from slider_isaac_lab.utils.slider_cfg import SliderCfg
from typing import Tuple, List


def get_intersection_interval_bounds(
    interval_1: torch.Tensor,  # [..., 2]
    interval_2: torch.Tensor   # [..., 2]
    ) -> torch.Tensor:             # [..., 2]
    """
    Returns intersection of two intervals if valid 

    Returns:
        torch.Tensor [..., 2]
    """
    # Check for NaNs
    if torch.isnan(interval_1).any() or torch.isnan(interval_2).any():
        return torch.full_like(interval_1, float('nan'))

    lower = torch.maximum(interval_1[..., 0], interval_2[..., 0])
    upper = torch.minimum(interval_1[..., 1], interval_2[..., 1])

    # Invalid intersections: intervals don’t overlap
    no_overlap = torch.logical_or(
        interval_1[..., 0] > interval_2[..., 1],
        interval_2[..., 0] > interval_1[..., 1]
    )

    lower[no_overlap] = float('nan')
    upper[no_overlap] = float('nan')

    return torch.stack((lower, upper), dim=-1)

def get_union_interval_bounds(
    interval_1: torch.Tensor, # [..., 2]
    interval_2:torch.Tensor # [..., 2]
    )->torch.Tensor: # [..., 2]
    """
    Returns union of two intervals 
    
    Returns:
        torch.Tensor [..., 2]
    """

    final_interval = torch.ones((*interval_1.shape[:-1], 2), dtype=interval_1.dtype, device=interval_1.device) * float('inf')

    final_interval[..., 0] = torch.minimum(interval_1[..., 0], interval_2[..., 0])
    final_interval[..., 1] = torch.maximum(interval_1[..., 1], interval_2[..., 1])
    
    return final_interval

def breaking_distance_based_action_bounds(
    m:torch.Tensor, 
    dt:torch.Tensor, 
    b_p:torch.Tensor,
    alpha:torch.Tensor, 
    v:torch.Tensor, 
    p:torch.Tensor, 
    f_max:torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:        
    """
    Defines analytical formulas for computing candidate control 
    intervals that satisfy state constraints.
    
    The state constraints are expressed as braking-distance limits
    relative to the outer boundary of a bounding box in the x-y plane of the world frame.    
    """
    
    # Candidate solutions from quadratic formula
    def get_root_1_minus(f_max, m, dt, b_p, gamma1, p, v):
        return -f_max*m*(-1/2*dt**2/m + dt*v/f_max \
            - dt*torch.sqrt(m*(2*f_max*(b_p*gamma1 + dt*v - gamma1*torch.abs(p) + p + torch.abs(p)) \
                - gamma1*m*v**2) + (1/4)*(f_max*dt - 2*m*v)**2)/(f_max*m))/dt**2

    def get_root_2_minus(f_max, m, dt, b_p, gamma1, p, v):
        return -f_max*m*(-1/2*dt**2/m + dt*v/f_max \
            + dt*torch.sqrt(m*(2*f_max*(b_p*gamma1 + dt*v - gamma1*torch.abs(p) + p + torch.abs(p))\
                - gamma1*m*v**2) + (1/4)*(f_max*dt - 2*m*v)**2)/(f_max*m))/dt**2

    def get_root_1_plus(f_max, m, dt, b_p, gamma1, p, v):
        return -f_max*m*((1/2)*dt**2/m + dt*v/f_max \
            - dt*torch.sqrt(-m*(2*f_max*(-b_p*gamma1 + dt*v + gamma1*torch.abs(p) + p - torch.abs(p))\
                + gamma1*m*v**2) + (1/4)*(f_max*dt + 2*m*v)**2)/(f_max*m))/dt**2

    def get_root_2_plus(f_max, m, dt, b_p, gamma1, p, v):
        return -f_max*m*((1/2)*dt**2/m + dt*v/f_max \
            + dt*torch.sqrt(-m*(2*f_max*(-b_p*gamma1 + dt*v + gamma1*torch.abs(p) + p - torch.abs(p))\
                + gamma1*m*v**2) + (1/4)*(f_max*dt + 2*m*v)**2)/(f_max*m))/dt**2
    
    r_m1 = get_root_1_minus(f_max, m, dt, b_p, alpha, p, v)
    r_m2 = get_root_2_minus(f_max, m, dt, b_p, alpha, p, v)
    lo_m = torch.minimum(r_m1, r_m2)
    hi_m = torch.maximum(r_m1, r_m2)

    r_p1 = get_root_1_plus(f_max, m, dt, b_p, alpha, p, v)
    r_p2 = get_root_2_plus(f_max, m, dt, b_p, alpha, p, v)
    lo_p = torch.minimum(r_p1, r_p2)
    hi_p = torch.maximum(r_p1, r_p2)
    
    interval_break_m = torch.stack([lo_m, hi_m], dim=-1)
    interval_break_p = torch.stack([lo_p, hi_p], dim=-1)

    return interval_break_m, interval_break_p

def get_abs_arg(m:torch.Tensor, dt:torch.Tensor,
                p:torch.Tensor, v:torch.Tensor
                )->torch.Tensor:
    return 2 * m * (-dt * v - p) / (dt ** 2)

def get_admissible_action_bounds(
    x: torch.Tensor,             
    x_bounds: torch.Tensor,      
    dt: float = 0.1,
    alpha: float = 0.05,) -> torch.Tensor:

    single = (x.ndim == 1)
    if single:
        x = x.unsqueeze(0)   

    device = x.device
    dtype = torch.float32

    x = x.to(dtype)
    x_bounds = x_bounds.to(device=device, dtype=dtype)

    m = torch.as_tensor(SliderCfg.mass, device=device, dtype=dtype)
    f_max = torch.as_tensor(SliderCfg.force_max, device=device, dtype=dtype)
    dt = torch.as_tensor(dt, device=device, dtype=dtype)
    alpha = torch.as_tensor(alpha, device=device, dtype=dtype)
    
    theta = x[:, 2]
    
    R = torch.stack((
        torch.stack((torch.cos(theta), -torch.sin(theta)), dim=-1),
        torch.stack((torch.sin(theta),  torch.cos(theta)), dim=-1)
        ), dim=-2).to(device=device, dtype=dtype)                         


    p_xy = x[:, 0:2]         
    v_xy = x[:, 3:5]          

    out = torch.full((x.shape[0], 2, 2), float('nan'), device=device, dtype=dtype)

    for i in range(2): 
        p = p_xy[:, i]
        v = v_xy[:, i]
        b_hi = x_bounds[i, 1]


        interval_break_m, interval_break_p = breaking_distance_based_action_bounds(m=m,dt=dt,b_p=b_hi,
                                                                                    alpha=alpha,
                                                                                    v=v,p=p,
                                                                                    f_max=f_max)
        arg_bound = get_abs_arg(m=m, 
                                dt=dt, 
                                p=p, 
                                v=v)

        arg_interval_m = torch.stack([torch.full_like(arg_bound, float('-inf')), arg_bound], dim=-1)
        arg_interval_p = torch.stack([arg_bound, torch.full_like(arg_bound, float('inf'))], dim=-1)
        
        # Actuator bounds are defined in local frame, convert to global 
        force_limit_interval = torch.tensor([-f_max, f_max]).to(R)
        force_limit_interval = torch.einsum('nij,k->nki', R, force_limit_interval)
        

        # Get candidate intervals
        m_interval = get_intersection_interval_bounds(interval_break_m, arg_interval_m)  # [..., 2]
        p_interval = get_intersection_interval_bounds(arg_interval_p, interval_break_p)  # [..., 2]

        valid_m = ~torch.isnan(m_interval).any(dim=-1)  # [...]
        valid_p = ~torch.isnan(p_interval).any(dim=-1)  # [...]

        final_interval = torch.full_like(m_interval, float('nan'))  # [..., 2]

        # only m valid
        mask = valid_m & ~valid_p
        final_interval[mask] = m_interval[mask]

        # only p valid
        mask = ~valid_m & valid_p
        final_interval[mask] = p_interval[mask]

        # both valid -> union
        both = valid_m & valid_p
        if both.any():
            final_interval[both] = get_union_interval_bounds(p_interval[both], m_interval[both])
        
        # Get final interval for admissible control based on state and actuator constraints 
        final_interval = get_intersection_interval_bounds(final_interval, force_limit_interval[..., i])
        
        out[:, i, :] = final_interval
        # out[:, i, 1] = final_interval[:, 1]

    return out[0] if single else out
