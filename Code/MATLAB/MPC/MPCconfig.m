function MPC = MPCconfig(cfg)

    MPC.Np = 20;          % prediction horizon
    MPC.Nc = 5;           % control horizon

    Mx_max = 10;
    My_max = 10;
    Mz_max = 10;

    MPC.u_min = cfg.u_min;
    MPC.u_max = cfg.u_max;

    %Weights
    MPC.Q = diag([10, 10, 50,  1, 1, 5]);
    MPC.R = diag([0.1, 0.1, 0.1]);

    cfg.x_op = [0 0 0 0 0 0];
    cfg.u_op = [0 0 0];

    [Al, Bl] = sliderLinearModel(cfg);

    Ts = cfg.Ts;
    [Ad,Bd] = discretizeBackwardEuler(Al,Bl,Ts);

    MPC.Ad = Ad;
    MPC.Bd = Bd;

    MPC.a_max = [ Mx_max;  My_max;  Mz_max];
    MPC.a_min = [-Mx_max; -My_max; -Mz_max];
end
function [A,B] = sliderLinearModel(cfg)

    % --- unpack params ---
    m   = cfg.mass;
    Izz = cfg.I_zz;

    x_op = cfg.x_op(:);   % ensure column

    theta0 = x_op(3);
    vx0    = x_op(4);
    vy0    = x_op(5);
    r0     = x_op(6);

    c = cos(theta0);
    s = sin(theta0);

    % --- A matrix (Jacobian df/dx at operating point) ---
    A = zeros(6,6);

    % row 1: x_dot = vx*c - vy*s
    A(1,3) = -vy0*c - vx0*s;   % d/dtheta
    A(1,4) =  c;               % d/dvx
    A(1,5) = -s;               % d/dvy

    % row 2: y_dot = vx*s + vy*c
    A(2,3) =  vx0*c - vy0*s;   % d/dtheta
    A(2,4) =  s;               % d/dvx
    A(2,5) =  c;               % d/dvy

    % row 3: theta_dot = r
    A(3,6) = 1;

    % row 4: vx_dot = r*vy + FxB/m
    A(4,5) =  r0;              % d/dvy
    A(4,6) =  vy0;             % d/dr

    % row 5: vy_dot = -r*vx + FyB/m
    A(5,4) = -r0;              % d/dvx
    A(5,6) = -vx0;             % d/dr


    % --- B matrix (df/du) ---
    B = zeros(6,3);
    B(4,1) = 1/m; 
    B(5,2) = 1/m; 
    B(6,3) = 1/Izz;
end
function [Ad,Bd] = discretizeBackwardEuler(A,B,Ts)
    I = eye(size(A));
    M = (I - Ts*A);
    Ad = M \ I;       % (I - Ts*A)^(-1)
    Bd = M \ (Ts*B);  % (I - Ts*A)^(-1) * Ts * B
end

