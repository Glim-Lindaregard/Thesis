clear all; close all; sldiagviewer.clearAll;


% --- Base config (immutable) ---
cfg0  = config();
MPC = MPCconfig(cfg0);
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
failureTime = [inf 2 inf 3 inf 1 inf 0.4];    %[1,2,3,4,5,6,7,8]

simTime = 15;

ref = [1,5,pi/2,0,0,0]';

init = [1 0 0 0 0 0];

% Constants
mass    = 4.436; % [kg]
I_zz = 1.092; % [kgm^2]


% PID parameters (TEMP)
Kp = diag([2 2 0.48]);
Kd = diag([8 8 1.2]);
ki = diag([0.05 0.05 0.064]);


% --- run ---
simout = sim('sliderSim');
%ad = [1,1,1.3]';

%[uOut,index,x] = findUfromAd_DA(ad,uCache(256).U,uCache(256).A);

%ap = uCache(256).A * uOut;
%Visualizion
if 0
    AniOpts.fps = 90;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 20;
    AniOpts.broken_len = 0.6;

    %animateTrajectory(simout,ref,cfg0,AniOpts,failureTime);
    
    %plotStates(simout,ref);
    
    %plotOtherStuff(simout);

    %---Visualize AMS facets---
AMSOpts = struct(...
    'FaceColor',        [0.8 0.95 0.95], ...
    'EdgeColor',        [0.3 0.3 0.3], ...
    'FaceAlpha',        0, ...
    'EdgeAlpha',        0, ...
    'LineWidth',        0.8, ...
    'BackgroundColor',  [1 1 1]*1, ...
    'GridColor',        0.3*[1 1 1], ...
    'GridAlpha',        1, ...
    'GridLineStyle',    '--', ...
    'FontSize',         30, ...
    'View',             [55 30], ...            % fixed camera (no spin)
    'Lighting',         false, ...              % off for consistent panels
    'ShowNormals',      false, ...
    'Fancy',            false, ...              % no orbit/camera spin
    'ShowProduced',     true, ...
    'ShowDesired',      true, ...
    'MaxPanels',        9, ...                  % up to 9 scenarios (8 thrusters)
    'TileOverride',     [], ...                 % e.g., [rows cols]; [] = auto
    'TitlePrefix',      'Scenario ', ...
    'ShowBasis',        true,...
    'Index',  index);


    %a = simout.logsout.get('a').Values.Data;
    %ad = simout.logsout.get('ad').Values.Data;


    %visualizeAMSGrid(uCache,failureTime,AMSOpts);



    visualizeAMS(uCache(256).U,uCache(256).A,1,AMSOpts,ad,ap,x)

    %visualizeSlider(cfg0,U);
end