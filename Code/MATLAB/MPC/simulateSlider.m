function [xHist, adHist, uHist, t] = simulateSlider(cfg,uCashe, failureTime)
% SIMULATE_SLIDER_MPC_ALLOC
%   Closed-loop sim: wrench-MPC + AMS-based allocation + failures.
%
% Inputs:
%   cfg            : struct with physics + geometry (used to build phys)
%   AMS_all        : cell array of AMS data per failure healthyIndex
%   A_all          : cell array of A matrices per failure healthyIndex
%   failurehealthyIndexFun : function handle @(k_or_t) -> healthyIndex index
%
% Outputs:
%   x_hist : 6xNt   state trajectory
%   ad_hist: 3x(Nt-1) desired wrench from MPC
%   u_hist : 8x(Nt-1) thruster commands from allocation
%   t      : 1xNt time vector

    % --- simulation parameters ---
    Ts       = cfg.Ts;
    simTime = cfg.simTime;
    x0       = cfg.x0;         % 6x1
    xRef    = cfg.xRef;      % 6x1

    t = 0:Ts:simTime;
    Nt = numel(t);

    % --- physics for MPC ---
    phys.m      = cfg.m;
    phys.Izz    = cfg.Izz;
    phys.FxMax = cfg.FxMax;   % rough global bounds
    phys.FyMax = cfg.FyMax;
    phys.TauMax = cfg.TauMax;

    % --- weights ---
    Q = cfg.Q;   % 6x6
    R = cfg.R;   % 3x3

    % --- build MPC ---
    mpc = mpc_init(Ts, cfg.N, Q, R, phys);

    % --- storage ---
    xHist  = zeros(6, Nt);
    xHist(:,1) = x0;

    adHist = zeros(3, Nt-1);
    uHist  = zeros(8, Nt-1);   % assuming 8 thrusters total

    % --- main loop ---
    for k = 1:Nt-1
        xk = xHist(:,k);

        % 1) MPC: current state -> desired wrench a_d(k)
        ad_k = mpcStep(mpc, xk, xRef);    % 3x1
        adHist(:,k) = ad_k;

        % 2) Failure healthyIndex at this step
        %    You decide whether failurehealthyIndexFun takes k or t(k)
        tk = t(k);
        healthyIndex = failureIndex(tk,failureTime);

        % 3) Get current A and AMS for this failure case
        ACurr   = uCashe(healthyIndex).A;     % 3x8 or 3x(#healthy)
        AMSCurr = uCashe(healthyIndex).U;   % whatever findUfromAd expects

        % 4) Allocation: desired wrench -> thruster commands
        %    Replace this line with your actual findUfromAd signature
        u_k = findUfromAd_DA(ad_k, AMSCurr,ACurr);
        % Make sure u_k comes back as 8x1 (pad zeros for failed thrusters if needed)
        uHist(:,k) = u_k;

        % 5) Real wrench applied to plant
        aReal = ACurr * u_k;       % 3x1 [Fx; Fy; tau] realized

        % 6) Plant update using the same discrete model as MPC
        xNext = full(mpc.F_RK4(xk, aReal));
        xHist(:,k+1) = xNext;
    end
end
function healthyIdx = failureIndex(t, failureTime)
% HEALTHYINDEXFROMFAILURETIMES
%   Convert per-thruster failure times -> healthy index for AMS/A_all.
%
% Inputs:
%   t           : scalar time [s]
%   failureTime : 1x8 vector, failureTime(i) = time when thruster i fails.
%                 Convention: t < failureTime(i) => healthy.
%                 Use Inf for "never fails".
%
% Output:
%   healthyIdx  : integer in [1..256].
%                 healthyIdx - 1 is the 8-bit healthyIndex
%                 (1=healthy, 0=failed, LSB=thruster 1).

    % logical 1x8: true if thruster i is healthy at time t
    isHealthy = (t < failureTime);    % t < t_fail => still healthy

    % convert boolean vector -> 8-bit healthyIndex with LSB = thruster 1
    healthyIndex = uint8(0);
    for thr = 1:8
        if isHealthy(thr)
            healthyIndex = bitor(healthyIndex, bitshift(uint8(1), thr-1));
        end
    end

    % convert to 1-based index
    healthyIdx = double(healthyIndex) + 1;
end

