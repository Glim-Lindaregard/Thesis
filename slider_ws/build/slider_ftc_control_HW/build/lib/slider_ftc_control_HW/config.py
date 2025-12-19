from dataclasses import dataclass
import numpy as np

### --- This file sets up a lot of parameters used throughout the package --- ###

# --- Dataclasses ---

@dataclass(frozen=True)
class MPCParams:
    Fx_max: float
    Fy_max: float
    Tau_max: float
    N: int
    Q: np.ndarray
    R: np.ndarray
    Ts: float


@dataclass(frozen=True)
class PhysicalParams:
    m: float
    Izz: float
    N_thrusters: int
    max_failed_thr: int
    u_min: np.ndarray
    u_max: np.ndarray
    A: np.ndarray


@dataclass(frozen=True)
class WorldLimits:
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    wall_buffer: float


@dataclass(frozen=True)
class SliderConfig:
    phys: PhysicalParams
    mpc: MPCParams
    world: WorldLimits



def make_default_config() -> SliderConfig:


    # --- Physical parameters ---
    m = 4.436
    Izz = 1.092
    N_thrusters = 8
    max_failed_thr = 4 #Changes how many AMS's are pre computed
    slider_radius = 0.2
    slider_table_length = 4
    max_thrust = 0.7


    A = np.array([
             [-1,     0,     1,     0,     1,     0,    -1,     0   ],
             [ 0,     1,     0,     1,     0,    -1,     0,    -1   ],
             [-0.14,  0.14,  0.14, -0.14, -0.14,  0.14,  0.14, -0.14],
         ], dtype=np.float64)


    # --- MPC parameters ---
    Ts = 0.1
    FxMax  = 1.5 * (2 * max_thrust)
    FyMax  = 1.5 * (2 * max_thrust)
    TauMax = 1.5 * (4 * max_thrust * slider_radius)

    N = 15
    Q = np.diag([8, 8, 4,4,8,0, 8, 8, 8])
    R = np.diag([0.0, 0.0, 0.0])

    mpc = MPCParams(
        Fx_max=FxMax,
        Fy_max=FyMax,
        Tau_max=TauMax,
        N=N,
        Q=Q,
        R=R,
        Ts=Ts
    )

    # thrust limits
    u_max = max_thrust * np.ones(N_thrusters, dtype=float)
    u_min = np.zeros(N_thrusters, dtype=float)


    phys = PhysicalParams(
        m=m,
        Izz=Izz,
        N_thrusters=N_thrusters,
        max_failed_thr = max_failed_thr,
        u_min=u_min,
        u_max=u_max,
        A=A,
    )

    # --- World limits ---
    wall_buffer = 2 * slider_radius
    x_min = -slider_table_length / 2 + wall_buffer
    x_max = slider_table_length / 2 - wall_buffer
    y_min = -slider_table_length / 2 + wall_buffer
    y_max = slider_table_length / 2 - wall_buffer

    world = WorldLimits(
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        wall_buffer=wall_buffer,
    )


    return SliderConfig(
        phys=phys,
        mpc=mpc,
        world=world,
    )
