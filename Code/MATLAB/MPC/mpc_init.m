function mpc = mpc_init(Ts, N, Q, R, phys)
% MPC_INIT_WRENCH  Build NMPC in wrench space a_d = [Fx; Fy; tau].
%
% Inputs:
%   Ts   : sampling time [s]
%   N    : prediction horizon
%   Q    : 6x6 state weight matrix
%   R    : 3x3 input weight matrix
%   phys : struct with fields m, Izz
%
% Output:
%   mpc  : struct with solver, bounds, dims, and F_RK4

    import casadi.*

    % --- dimensions ---
    nStates   = 6;
    nControls = 3;

    % --- symbolic states and input ---
    x     = SX.sym('x');
    y     = SX.sym('y');
    theta = SX.sym('theta');
    vx    = SX.sym('vx');
    vy    = SX.sym('vy');
    r     = SX.sym('r');

    states  = [x; y; theta; vx; vy; r];

    ad     = SX.sym('ad', nControls, 1);   % [Fx; Fy; tau]

    % --- parameters ---
    m   = phys.m;
    Izz = phys.Izz;

    Fx  = ad(1);
    Fy  = ad(2);
    Tau = ad(3);

    % world-frame accel
    fxx = Fx*cos(theta) - Fy*sin(theta);
    fyy = Fx*sin(theta) + Fy*cos(theta);


    dwdt = [vx;
            vy;
            r;
            fxx/m;
            fyy/m;
            Tau/Izz];

    f = Function('f', {states, ad}, {dwdt});

    % --- RK4 discrete-time model ---
    XkSym    = SX.sym('Xk', nStates);
    UkSym    = SX.sym('Uk', nControls);
    XnextSym = RK4_step(f, XkSym, UkSym, Ts);

    F_RK4 = Function('F_RK4', {XkSym, UkSym}, {XnextSym});

    % --- optimization variables ---
    U = SX.sym('U', nControls, N);      % 3xN
    X = SX.sym('X', nStates, N+1);      % 6x(N+1)

    % Parameters P = [x0; x_ref]
    P = SX.sym('P', 2*nStates, 1);

    % --- objective & constraints ---
    obj = 0;
    g   = [];

    % initial condition: X(:,1) = x0
    g = [g; X(:,1) - P(1:nStates)];

    xRefSym = P(nStates+1:end);

    for k = 1:N
        st   = X(:,k);
        con  = U(:,k);
        st_n = X(:,k+1);

        % dynamics constraint
        st_n_pred = F_RK4(st, con);
        g = [g; st_n - st_n_pred];

        % tracking error
        err = st - xRefSym;

        % cost
        obj = obj + err.'*Q*err + con.'*R*con;
    end

    % stack decision variables
    OPTvariables = [reshape(X, nStates*(N+1), 1);
                     reshape(U, nControls*N,     1)];

    nlp_prob = struct('f', obj, 'x', OPTvariables, 'g', g, 'p', P);

    % --- solver options ---
    opts = struct;
    opts.ipopt.max_iter = 100;
    opts.ipopt.print_level = 0;
    opts.print_time = false;
    opts.ipopt.warm_start_init_point = 'no';

    solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

    % --- bounds ---
    lbg = zeros(length(g),1);
    ubg = zeros(length(g),1);

    lbx = [];
    ubx = [];

    xMin = phys.xMin;
    yMin = phys.yMin;
    xMax = phys.xMax;
    yMax = phys.yMax;


    %state bounds (unconstrained)

    for i = 1:(N+1)
        lbx = [lbx; -inf*ones(nStates,1)];
        ubx = [ubx;  inf*ones(nStates,1)];
        idx = stateIndex(i);
        lbx(idx.x) = xMin;
        ubx(idx.x) = xMax;
        lbx(idx.y) = yMin;
        ubx(idx.y) = yMax;
    end

    
    % control bounds (rough symmetric box for Fx,Fy,Tau)
    FxMax  = phys.FxMax;
    FyMax  = phys.FyMax;
    TauMax = phys.TauMax;
    
    for i = 1:N
        lbx = [lbx; -[FxMax; FyMax; TauMax]];
        ubx = [ubx;  [FxMax; FyMax; TauMax]];
    end

    % --- pack outputs ---
    mpc.solver       = solver;
    mpc.lbx          = lbx;
    mpc.ubx          = ubx;
    mpc.lbg          = lbg;
    mpc.ubg          = ubg;
    mpc.OPTvariablesLen = numel(OPTvariables);

    mpc.nStates   = nStates;
    mpc.nControls = nControls;
    mpc.N          = N;
    mpc.Ts         = Ts;

    mpc.F_RK4      = F_RK4;
    mpc.Q          = Q;
    mpc.R          = R;

end

% --- local RK4 helper ---
function X_next = RK4_step(f, Xk, Uk, Ts)
    k1 = f(Xk,          Uk);
    k2 = f(Xk + Ts/2*k1, Uk);
    k3 = f(Xk + Ts/2*k2, Uk);
    k4 = f(Xk + Ts   *k3, Uk);
    X_next = Xk + Ts/6*(k1 + 2*k2 + 2*k3 + k4);
end
function idx = stateIndex(i)
    nrStates = 6; base = (i-1)*nrStates;
    idx.x     = base +1;
    idx.y     = base +2;
    idx.theta = base + 3;
    idx.vx    = base + 4;
    idx.vy    = base + 5;
    idx.tau   = base + 6;
end