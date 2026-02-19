clc; clear; close all;
addpath('~/workspace/MATLAB_Packages');

%Set when which thruster fails
failureTime = [inf inf inf inf inf inf inf inf];    %[1,2,3,4,5,6,7,8]

%Build and store all AMSs and A's
cfg0 = config();
uCache = buildAllAMS(cfg0);  

% --- run simulation ---


%[xHist, adHist, uHist,aRealHist, t] = simulateSlider(cfg0,uCache,failureTime);
%Visualizion
if 0

    animateTrajectory(t,xHist,uHist,cfg0,failureTime);
    
    plotStates(t,xHist,cfg0.xRef);
    %plotCopare(t,xHist,xHistC,cfg0.xRef);
    
    %plotOtherStuff(t,adHist,aRealHist);

    %visualizeAMSGrid(uCache,failureTime,AMSOpts);

    %visualizeAMS(uCache(256).U,uCache(256).A,'normal');

    %visualizeSlider(cfg0,U);
end

function uCache = buildAllAMS(cfg0)
    % --- Base config ---
    N     = cfg0.N_thrusters;                 % should be 8
     
    %Cashe to store all 256 AMSs and As
    uCache = struct([]);
    
    %Build all 256 AMSs
    for k = 0:(2^N - 1)
        mask = logical(bitget(k, 1:N));   % 1 = healthy, 0 = failed
        
        cfgi = cfg0;
        cfgi.A(:, ~mask) = 0;             % zero failed thruster columns
        cfgi.u_min(~mask) = 0;            % zero min/max ranges
        cfgi.u_max(~mask) = 0;
        cfgi.N = N;                       % keep same dimension
        
        Ui= buildAMS_row(cfgi);
        uCache(k+1).U = Ui;
        uCache(k+1).A = cfgi.A;
        uCache(k+1).mask = mask;
    end
end