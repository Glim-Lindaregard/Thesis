clc; clear; close all;

% --- Base config (immutable) ---
cfg0  = config();
m     = size(cfg0.A,2);                  % should be 8

uCache = struct([]);

global ads;

ads = [];


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
failureTime = [1 inf inf inf inf inf 2 inf];    %[1,2,3,4,5,6,7,8]


% --- run ---
%simout = sim('sliderSim');
[xHist, adHist, uHist, t] = simulateSlider(cfg0,uCache,failureTime);

%Visualizion
if 1
    AniOpts.fps = 90;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 20;
    AniOpts.broken_len = 0.6;

    animateTrajectory(t,xHist,uHist,cfg0,AniOpts,failureTime);
    
    plotStates(t,xHist,cfg0.xRef);
    
    %plotOtherStuff(simout);

    %---Visualize AMS facets---
    AMSOptsExport = struct(...
        'FaceColor',        [0.8 0.95 0.95], ...
        'EdgeColor',        [0.3 0.3 0.3], ...
        'FaceAlpha',        0.95, ...
        'EdgeAlpha',        1, ...
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
        'Index',  1);

    AMSOpts = struct(...
        'FaceColor',        [0.8 0.95 0.95], ...
        'EdgeColor',        [0.3 0.3 0.3]*3, ...
        'FaceAlpha',        0.7, ...
        'EdgeAlpha',        0.2, ...
        'LineWidth',        3, ...
        'BackgroundColor',  [1 1 1]*0.2, ...
        'GridColor',        0.2*[1 1 1], ...
        'GridAlpha',        1, ...
        'GridLineStyle',    '--', ...
        'FontSize',         11, ...
        'View',             [55 30], ...            % fixed camera (no spin)
        'Lighting',         false, ...              % off for consistent panels
        'ShowNormals',      false, ...
        'Fancy',            false, ...              % no orbit/camera spin
        'ShowProduced',     false, ...
        'ShowDesired',      true, ...
        'MaxPanels',        9, ...                  % up to 9 scenarios (8 thrusters)
        'TileOverride',     [], ...                 % e.g., [rows cols]; [] = auto
        'TitlePrefix',      'Scenario ', ...
        'ShowBasis',        false,...
        'Index',  1);


    %a = simout.logsout.get('a').Values.Data;
    %ad = simout.logsout.get('ad').Values.Data;


    %visualizeAMSGrid(uCache,failureTime,AMSOpts);


    %visualizeAMS(uCache(256).U,uCache(256).A,0,AMSOpts,0,ads,0);


    %visualizeSlider(cfg0,U);
end