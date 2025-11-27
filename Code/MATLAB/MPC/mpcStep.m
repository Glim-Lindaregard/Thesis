function ad_k = mpcStep(mpc, xCurrent, xRef, uBounds)
% MPC_STEP  One MPC step: current state -> desired wrench a_d,
%           with time-varying box constraints from AMS.
%
% Inputs:
%   mpc      : struct from mpc_init
%   xCurrent : 6x1 current state
%   xRef     : 6x1 reference state
%   uBounds  : struct with fields
%              Fx_min, Fx_max, Fy_min, Fy_max, Tau_min, Tau_max
%
% Output:
%   ad_k     : 3x1 desired wrench [Fx; Fy; tau] at current step

    import casadi.*

    nStates   = mpc.nStates;
    nControls = mpc.nControls;   % should be 3
    N         = mpc.N;

    % pack parameters P = [x0; x_ref]
    P = [xCurrent; xRef];

    % initial guess (cold start for now)
    x0opt = zeros(mpc.OPTvariablesLen, 1);

    % --- build per-step bounds for controls from uBounds ---
    lbx = mpc.lbx;   % start from template
    ubx = mpc.ubx;

    % Index where controls start in OPTvariables:
    % [ X(0..N) ; U(0..N-1) ]
    nX = nStates*(N+1);

    lb_u_stage = [uBounds.Fx_min;  uBounds.Fy_min;  uBounds.Tau_min];
    ub_u_stage = [uBounds.Fx_max;  uBounds.Fy_max;  uBounds.Tau_max];

    % Apply same box for all N control stages
    for i = 1:N
        idx = nX + (i-1)*nControls + (1:nControls);
        lbx(idx) = lb_u_stage;
        ubx(idx) = ub_u_stage;
    end

    % --- solve NLP with updated bounds ---
    sol = mpc.solver('x0', x0opt, ...
                     'lbx', lbx, 'ubx', ubx, ...
                     'lbg', mpc.lbg, 'ubg', mpc.ubg, ...
                     'p',  P);

    solX = full(sol.x);

    % extract first control input a_d(k)
    U_sol = reshape(solX(nX+1:end), nControls, N);
    ad_k  = U_sol(:,1);
end
