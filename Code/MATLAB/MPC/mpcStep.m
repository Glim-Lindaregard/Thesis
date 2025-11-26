function ad_k = mpcStep(mpc, xCurrent, xRef)
% MPC_STEP_WRENCH  One MPC step: current state -> desired wrench a_d.
%
% Inputs:
%   mpc       : struct from mpc_init_wrench
%   x_current : 6x1 current state
%   x_ref     : 6x1 reference state
%
% Output:
%   ad_k      : 3x1 desired wrench [Fx; Fy; tau] at current step

    import casadi.*

    nStates   = mpc.nStates;
    nControls = mpc.nControls;
    N          = mpc.N;

    % pack parameters P = [x0; x_ref]
    P = [xCurrent; xRef];

    % initial guess (cold start for now)
    x0opt = zeros(mpc.OPTvariablesLen, 1);

    % solve NLP
    sol = mpc.solver('x0', x0opt, ...
                     'lbx', mpc.lbx, 'ubx', mpc.ubx, ...
                     'lbg', mpc.lbg, 'ubg', mpc.ubg, ...
                     'p',  P);

    solX = full(sol.x);

    % extract first control input a_d(k)
    nX    = nStates*(N+1);
    U_sol = reshape(solX(nX+1:end), nControls, N);
    ad_k  = U_sol(:,1);
end
