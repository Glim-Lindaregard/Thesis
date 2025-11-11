clear all; close all;

% --- Base config (immutable) ---
cfg0  = config();
m     = size(cfg0.A,2);                  % should be 8

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



%healthy = bin2dec('00111111') + 1;  % [8,7,6,5,4,3,2,1]
failureTime = [inf 2 inf inf inf 1 inf inf];    %[1,2,3,4,5,6,7,8]

simTime = 15;

ref = [1,1.3,0.1,0,0,0]';

init = [1 0.1 0.1 0 0 0];

% Constants
mass    = 4.436; % [kg]
I_zz = 1.092; % [kgm^2]


% PID parameters (TEMP)
Kp = diag([2 2 0.48]);
Kd = diag([8 8 1.2]);
ki = diag([0.05 0.05 0.064]);


% --- run ---
simout = sim('sliderSim');


%Visualizion
if 1
    AniOpts.fps = 90;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 20;
    AniOpts.broken_len = 0.6;

    animateTrajectory(simout,ref,cfg0,AniOpts,failureTime);
    
    plotStates(simout,ref);
    
    plotOtherStuff(simout);

    %---Visualize AMS facets---
AMSOpts = struct(...
    'FaceColor',        [0.78 0.92 0.92], ...
    'EdgeColor',        [0.15 0.15 0.18], ...
    'FaceAlpha',        0.95, ...
    'EdgeAlpha',        0.85, ...
    'LineWidth',        0.7, ...
    'BackgroundColor',  [1 1 1]*0.2, ...
    'GridColor',        [0.65 0.65 0.70], ...
    'GridAlpha',        0.35, ...
    'GridLineStyle',    '--', ...
    'FontSize',         11, ...
    'View',             [45 30], ...            % fixed camera (no spin)
    'Lighting',         false, ...              % off for consistent panels
    'ShowNormals',      false, ...
    'Fancy',            false, ...              % no orbit/camera spin
    'ShowProduced',     true, ...
    'ShowDesired',      true, ...
    'MaxPanels',        9, ...                  % up to 9 scenarios (8 thrusters)
    'TileOverride',     [], ...                 % e.g., [rows cols]; [] = auto
    'TitlePrefix',      'Scenario ');


    a = simout.logsout.get('a').Values.Data;
    ad = simout.logsout.get('ad').Values.Data;


    visualizeAMSGrid(uCache,failureTime,AMSOpts);

    %visualizeSlider(cfg,U);
end