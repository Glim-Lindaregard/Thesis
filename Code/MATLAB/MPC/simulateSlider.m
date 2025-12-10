function [xHist, adHist, uHist,aRealHist, t] = simulateSlider(cfg,uCashe, failureTime)
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
    phys.xMin = cfg.xMin;
    phys.xMax = cfg.xMax;
    phys.yMin = cfg.yMin;
    phys.yMax = cfg.yMax;

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

    aRealHist = zeros(3, Nt-1);

    % --- main loop ---
    for k = 1:Nt-1
        xk = xHist(:,k);
        tk = t(k);
    
        % --- 1) Failure pattern & current A / AMS ---
        healthyIndex = failureIndex(tk, failureTime);   % 1..256
    
        ACurr   = uCashe(healthyIndex).A;    % 3 x 8 (with zeroed failed cols)
        AMSCurr = uCashe(healthyIndex).U;    % AMS facets for this failure mode
        mask    = uCashe(healthyIndex).mask; % 1 x 8 logical, healthy thrusters
    
        % --- 2) Compute admissible wrench box from AMS / A ---
        % cfg.u_min / cfg.u_max are the *nominal* thruster limits
        uBounds = computeBounds(ACurr, cfg.u_min, cfg.u_max, mask);

        %xRef = cfg.refFun(tk);

        xRef = cfg.xRef;
    
        % --- 3) MPC: current state -> desired wrench a_d(k) within bounds ---
        ad_k = mpcStep(mpc, xk, xRef, uBounds);    % 3x1
        adHist(:,k) = ad_k;
    
        % --- 4) Allocation: desired wrench -> thruster commands (DA) ---
        u_k = findUfromAd_DA(ad_k, AMSCurr, ACurr);
        uHist(:,k) = u_k;
    
        % --- 5) Real wrench applied to plant ---
        aReal = ACurr * u_k;       % 3x1 [Fx; Fy; tau] realized
        aRealHist(:,k) = aReal;    % LOG
    
        % --- 6) Plant update using same discrete model as MPC ---
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

