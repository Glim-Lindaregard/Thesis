close all;

% --- Base config (immutable) ---
cfg0  = config();
A0    = cfg0.A; 
umin0 = cfg0.u_min; 
umax0 = cfg0.u_max;
m     = size(A0,2);                  % should be 8
MASK_ALL = uint16(2^m - 1);          % 255 for m=8

% --- Cache (256 slots: mask 0..255 maps to 1..256) ---
AMS_cache = zeros(256,1);

% Healthy
AMS_H = buildAMS_row(cfg0);
AMS_H.meta.mask        = MASK_ALL;
AMS_H.meta.active_idx  = 1:m;
AMS_cache(double(MASK_ALL)+1) = AMS_H;

% Single-failures
for i = 1:m
    idx   = setdiff(1:m,i);          % active thrusters
    cfgi  = cfg0;
    cfgi.A     = A0(:,idx);
    cfgi.u_min = umin0(idx);
    cfgi.u_max = umax0(idx);
    cfgi.N = size(cfgi.A,2);

    AMSi = buildAMS_row(cfgi);
    mask_i = bitset(MASK_ALL, i, 0);

    AMSi.meta.mask       = mask_i;
    AMSi.meta.active_idx = idx;
    AMS_cache(double(mask_i)+1) = AMSi;
end

% --- Using it at runtime ---
% health: 1=healthy, 0=failed (LSB=thruster 1)
health = [1 1 1 1 1 1 1 1];                               % example
mask = uint16(sum(uint16(health).*uint16(2.^(0:m-1))));
AMS = AMS_cache(double(mask)+1);          % fetch




simTime = 10;

ref = [1,2,pi/2,0,0,0]';

init = [1 0 pi/2 0 0 0];

% Constants
m    = 4.436; % [kg]
I_zz = 1.092; % [kgm^2]


failureTime = 0.2;


% PID parameters (TEMP)
Kp = diag([2.5 2.5 0.48]);
Kd = diag([6 6 1.2]);
ki = diag([0.05 0.05 0.064]);


% --- run ---
simout = sim('sliderSim');


%Visualizion
if 1
    AniOpts.fps = 250;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 60;

    %animateTrajectory(simout,ref,AniOpts);
    
    %plotStates(simout,ref);
    
    %plotOtherStuff(simout);

    %---Visualize AMS facets---
    AMSopts = struct( ...
        'FaceColor', [0.78 1 0.92], ... % pastel cyan
        'EdgeColor', 1*[1 1 1], ...
        'FaceAlpha', 0.95, ...
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
        'ShowProduced', false,...
        'ShowDesired', false,...
        'ShowBasis', false);


    %visualizeAMS(AMS_current,AMSopts);

    %visualizeSlider(cfg,U);
end


%Helper functions
function result = isHealthy(mask,i)
    result = bitget(mask,i);
end

function mask = setBroken(mask,i)
    mask = bitset(mask,i,0);
end