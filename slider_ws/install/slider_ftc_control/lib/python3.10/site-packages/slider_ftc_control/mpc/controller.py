from dataclasses import dataclass
import numpy as np
import casadi as ca
from slider_ftc_control.config import SliderConfig
from .bounds import WrenchBounds
from tf_transformations import euler_from_quaternion


@dataclass
class MPCDims:
    n_states: int
    n_controls: int
    N: int
    Ts: float
    opt_len: int
    nX: int

# --- Controller based on Sathyas design --- 
class MPCController:
    def quat_conj(self,q):
    # q = [x,y,z,w]
        return ca.vertcat(-q[0], -q[1], -q[2], q[3])

    def rotate_vec_by_quat(self,q, v):
        vq = ca.vertcat(v[0], v[1], v[2], 0)
        return self.quat_mul(self.quat_mul(q, vq), self.quat_conj(q))[0:3]
    
    def __init__(self, cfg: SliderConfig):
        self.cfg = cfg
        self.prev_solution = None # Used as initial guess for solver

        Ts = cfg.mpc.Ts
        N = cfg.mpc.N
        Q = cfg.mpc.Q
        R = cfg.mpc.R


        m = cfg.phys.m
        Izz = cfg.phys.Izz

        xMin = cfg.world.x_min
        xMax = cfg.world.x_max
        yMin = cfg.world.y_min
        yMax = cfg.world.y_max

        # --- Default mpc bounds ---
        FxMax = cfg.mpc.Fx_max
        FyMax = cfg.mpc.Fy_max
        TauMax = cfg.mpc.Tau_max

        # --- dimensions ---
        n_states = 9
        n_controls = 3

        # --- symbolic states and inputs ---
        x = ca.SX.sym('x')
        y = ca.SX.sym('y')
        q = ca.SX.sym('q', 4) # quaternion
        #theta = ca.SX.sym('theta')
        vx = ca.SX.sym('vx')
        vy = ca.SX.sym('vy')
        r = ca.SX.sym('r')

        states = ca.vertcat(x, y, q, vx, vy, r)  # 7x1

        ad = ca.SX.sym('ad', n_controls, 1)  # [Fx; Fy; Tau]

        Fx = ad[0]
        Fy = ad[1]
        Tau = ad[2]
        f_body = ca.vertcat(Fx, Fy, 0)

        f_world = self.rotate_vec_by_quat(q, f_body)
        fxx = f_world[0]
        fyy = f_world[1]

        omega_q = ca.vertcat(0, 0, r, 0)
        q_dot = 0.5 * self.quat_mul(q, omega_q)

        dwdt = ca.vertcat(
            vx,
            vy,
            q_dot,  # dq/dt
            fxx / m,
            fyy / m,
            Tau / Izz,
        )

        f = ca.Function('f', [states, ad], [dwdt])

        # --- RK4 discrete model ---
        Xk_sym = ca.SX.sym('Xk', n_states)
        Uk_sym = ca.SX.sym('Uk', n_controls)

        def rk4_step(f, Xk, Uk, Ts):
            k1 = f(Xk,           Uk)
            k2 = f(Xk + Ts/2*k1, Uk)
            k3 = f(Xk + Ts/2*k2, Uk)
            k4 = f(Xk + Ts   *k3, Uk)
            return Xk + Ts/6 * (k1 + 2*k2 + 2*k3 + k4)

        Xnext_sym = rk4_step(f, Xk_sym, Uk_sym, Ts)
        F_RK4 = ca.Function('F_RK4', [Xk_sym, Uk_sym], [Xnext_sym])

        # --- optimization variables ---
        U = ca.SX.sym('U', n_controls, N)      # 3xN
        X = ca.SX.sym('X', n_states, N + 1)    # 6x(N+1)

        # --- Parameters P = [x0; x_ref] ---
        P = ca.SX.sym('P', 2 * n_states, 1)

        obj = 0
        g = []
        x0_param = P[0:n_states]
        xRef_sym = P[n_states:2*n_states]

        g.append(X[:, 0] - x0_param)

        # --- Build dynamics contraints (x_next - x_pred(RK4) = 0)and objective function xQx^T + uRu^T
        for k in range(N):
            st = X[:, k]
            con = U[:, k]
            st_n = X[:, k + 1]
            q=st[2:6]
            q_d=xRef_sym[2:6]
            q_err=self.quat_mul(q,self.quat_conj(q_d))


            st_n_pred = F_RK4(st, con)
            g.append(st_n - st_n_pred)

            err = st - xRef_sym
            err[2:6] = q_err  # Quaternion error
            obj = obj + ca.mtimes([err.T, Q, err]) + ca.mtimes([con.T, R, con])

        g_stack = ca.vertcat(*g)

        # --- One long (n_states*(N+1) + n_controls*N) x 1  opt vector ---
        OPTvariables = ca.vertcat(
            ca.reshape(X, n_states * (N + 1), 1),
            ca.reshape(U, n_controls * N, 1),
        )

        nlp_prob = {
            'f': obj,
            'x': OPTvariables,
            'g': g_stack,
            'p': P,
        }

        opts = {
            'ipopt': {
                'max_iter': 100,
                'print_level': 0,
                'warm_start_init_point': 'no',
            },
            'print_time': False,
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp_prob, opts)

        # --- bounds on g (g == 0) AKA X_{k+1} - RK4(X_k,U_k) = 0 ---
        lbg = np.zeros(g_stack.shape[0])
        ubg = np.zeros(g_stack.shape[0])

        # --- bounds on x ---
        lbx = []
        ubx = []

        # --- state bounds ---
        for i in range(N + 1):
            # default: unconstrained
            lbx.extend([-np.inf] * n_states)
            ubx.extend([+np.inf] * n_states)

            base = i * n_states
            idx_x = base + 0
            idx_y = base + 1
            # Table bounds
            lbx[idx_x] = xMin
            ubx[idx_x] = xMax
            lbx[idx_y] = yMin
            ubx[idx_y] = yMax

        # --- Attatch control bounds ---
        for _ in range(N):
            lbx.extend([ -FxMax,  -FyMax,  -TauMax])
            ubx.extend([ +FxMax,  +FyMax,  +TauMax])

        lbx = np.array(lbx, dtype=float)
        ubx = np.array(ubx, dtype=float)

        nX = n_states * (N + 1)
        opt_len = OPTvariables.shape[0]

        self.solver = solver
        self.lbg = lbg
        self.ubg = ubg
        self.lbx_template = lbx
        self.ubx_template = ubx
        self.F_RK4 = F_RK4

        self.dims = MPCDims(
            n_states=n_states,
            n_controls=n_controls,
            N=N,
            Ts=Ts,
            opt_len=int(opt_len),
            nX=nX,
        )

    # ------------------------------------------------------------------    
    def quat_mul(self,q1, q2):
        x1,y1,z1,w1 = q1[0],q1[1],q1[2],q1[3]
        x2,y2,z2,w2 = q2[0],q2[1],q2[2],q2[3]
        return ca.vertcat(
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
        w1*w2 - x1*x2 - y1*y2 - z1*z2
        )

    def step(self,
             x_current: np.ndarray,
             x_ref: np.ndarray,
             bounds: WrenchBounds) -> np.ndarray:

        n_states = self.dims.n_states
        n_controls = self.dims.n_controls
        N = self.dims.N
        nX = self.dims.nX

        x_current = np.asarray(x_current, float).reshape(n_states)
        x_ref = np.asarray(x_ref, float).reshape(n_states)

        # --- Parameters P = [x0; x_ref] ---
        P = np.concatenate([x_current, x_ref])

        # --- Initial guess (cold start) initial guess is 0,0,0,0,0,0 ---
        #x0opt = np.zeros(self.dims.opt_len)

        # --- Initial guess (warm start) initial guess is previous state ---
        x0opt = self.warm_start()

        lbx = self.lbx_template.copy()
        ubx = self.ubx_template.copy()

        # --- Attatch custom bounds from "bounds" depending on thruster health ---
        lb_u_stage = np.array(
            [bounds.Fx_min, bounds.Fy_min, bounds.Tau_min],
            dtype=float,
        )
        ub_u_stage = np.array(
            [bounds.Fx_max, bounds.Fy_max, bounds.Tau_max],
            dtype=float,
        )

        for i in range(N):
            idx = nX + i * n_controls
            lbx[idx:idx + n_controls] = lb_u_stage
            ubx[idx:idx + n_controls] = ub_u_stage

        sol = self.solver(
            x0=x0opt,
            lbx=lbx,
            ubx=ubx,
            lbg=self.lbg,
            ubg=self.ubg,
            p=P,
        )

        # --- Log solution to use for initial guess next step ---
        self.prev_solution = sol['x'].full().flatten()

        solX = np.array(sol['x']).reshape(-1)

        # --- Extract first control input a_d(k) ---
        U_sol = solX[nX:]
        U_mat = U_sol.reshape(N, n_controls).T  # shape (3, N)
        ad_k = U_mat[:, 0]

        return ad_k

    # --- Helper to construct initial guess from previous state ---
    def warm_start(self) -> np.ndarray:

        n_states   = self.dims.n_states
        n_controls = self.dims.n_controls
        N          = self.dims.N
        nX         = self.dims.nX
        opt_len    = self.dims.opt_len


        if self.prev_solution is None or self.prev_solution.shape[0] != opt_len:
            return np.zeros(opt_len, dtype=float)

        prev = self.prev_solution
        x0opt = np.zeros_like(prev)

        for i in range(N):
            x0opt[i * n_states : (i + 1) * n_states] = \
                prev[(i + 1) * n_states : (i + 2) * n_states]

        x0opt[N * n_states : (N + 1) * n_states] = \
            prev[N * n_states : (N + 1) * n_states]

        ctrl_start = nX

        for i in range(N - 1):
            x0opt[ctrl_start + i * n_controls : ctrl_start + (i + 1) * n_controls] = \
                prev[ctrl_start + (i + 1) * n_controls : ctrl_start + (i + 2) * n_controls]

        x0opt[ctrl_start + (N - 1) * n_controls : ctrl_start + N * n_controls] = \
            prev[ctrl_start + (N - 1) * n_controls : ctrl_start + N * n_controls]

        return x0opt

