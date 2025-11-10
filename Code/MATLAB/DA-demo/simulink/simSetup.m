clear all; close all;

% --- Base config (immutable) ---
cfg0  = config();
m     = size(cfg0.A,2);                  % should be 8

U0 = Copy_of_buildAMS_row(cfg0);

% uCashe(m+1).U = U0;
% uCashe(m+1).A = cfg0.A;

uCache = struct([]);


for k = 0:(2^m - 1)
    mask = logical(bitget(k, 1:m));   % 1 = healthy, 0 = failed
    
    cfgi = cfg0;
    cfgi.A(:, ~mask) = 0;             % zero failed thruster columns
    cfgi.u_min(~mask) = 0;            % zero min/max ranges
    cfgi.u_max(~mask) = 0;
    cfgi.N = m;                       % keep same dimension
    
    Ui= Copy_of_buildAMS_row(cfgi);
    uCache(k+1).U = Ui;
    uCache(k+1).A = cfgi.A;
    uCache(k+1).mask = mask;
end


healthy = bin2dec('10011111') + 1;
failureTime = 0;

simTime = 10;

ref = [3,5,0,0,0,0]';

init = [3 0 0 0 0 0];

% Constants
mass    = 4.436; % [kg]
I_zz = 1.092; % [kgm^2]





% PID parameters (TEMP)
Kp = diag([2 2 0.48]);
Kd = diag([8 8 1.2]);
ki = diag([0.05 0.05 0.064]);


% --- run ---
simout = sim('sliderSim');

m = size(cfg0.A,2);
%assert(ismember(broken,1:m+1),'broken must be 1..m+1');

% A_used   = uCache(broken).A;
% zeroCols = find(all(abs(A_used)<1e-12,1));   % which thruster(s) are zeroed
% disp(struct('broken',broken,'zeroed_columns',zeroCols))



%Visualizion
if 1
    AniOpts.fps = 500;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 200;

    animateTrajectory(simout,ref,cfg0,AniOpts,9,failureTime);
    
    plotStates(simout,ref);
    
    plotOtherStuff(simout);

    %---Visualize AMS facets---
    AMSopts = struct( ...
        'FaceColor', [0.78 1 0.92], ... % pastel cyan
        'EdgeColor', 1*[1 1 1], ...
        'FaceAlpha', 0.1, ...
        'EdgeAlpha', 0.8, ...
        'LineWidth', 0.6, ...
        'BackgroundColor', 0.2*[1 1 1], ...
        'GridColor', 0.4*[0.9 0.9 0.9], ...
        'GridAlpha', 1, ...
        'GridLineStyle', '-', ...
        'UseOctantColors', false, ...
        'Lighting', true, ...
        'ShowNormals', false,...
        'fps', 10, ...
        'Fancy', false, ...
        'Index', 0, ...
        'ShowProduced', true,...
        'ShowDesired', true,...
        'ShowBasis', false);
    a = simout.logsout.get('a').Values.Data;
    ad = simout.logsout.get('ad').Values.Data;




    visualizeAMS(uCache(healthy).U,uCache(healthy).A,1,AMSopts,a(:,:,10),ad(:,:,10),1);

    %visualizeSlider(cfg,U);
end