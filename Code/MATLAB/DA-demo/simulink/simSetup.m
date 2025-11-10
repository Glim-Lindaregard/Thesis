clear all; close all;

% --- Base config (immutable) ---
cfg0  = config();
m     = size(cfg0.A,2);                  % should be 8

U0 = Copy_of_buildAMS_row(cfg0);

uCashe(m+1).U = U0;
uCashe(m+1).A = cfg0.A;

% Single-failures
for i = 1:2^m
    cfgi  = cfg0;
    cfgi.A(:,i) = 0;
    cfgi.u_min(i) = 0;
    cfgi.u_max(i) = 0;
    cfgi.N = size(cfgi.A,2);

    Ui = Copy_of_buildAMS_row(cfgi);
    uCashe(i).U = Ui;
    uCashe(i).A = cfgi.A;
end


broken = 4;
failureTime = 1;

simTime = 20;

ref = [3,5,pi + 0.3,0,0,0]';

init = [0 0 pi/2 0 0 0];

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
assert(ismember(broken,1:m+1),'broken must be 1..m+1');

A_used   = uCashe(broken).A;
zeroCols = find(all(abs(A_used)<1e-12,1));   % which thruster(s) are zeroed
disp(struct('broken',broken,'zeroed_columns',zeroCols))



%Visualizion
if 1
    AniOpts.fps = 500;
    AniOpts.saveVideo = false;
    AniOpts.videoName = 'sliderAnimation.mp4';
    AniOpts.simSpeed = 200;

    animateTrajectory(simout,ref,cfg0,AniOpts,broken,failureTime);
    
    plotStates(simout,ref);
    
    %plotOtherStuff(simout);

    %---Visualize AMS facets---
    AMSopts = struct( ...
        'FaceColor', [0.78 1 0.92], ... % pastel cyan
        'EdgeColor', 1*[1 1 1], ...
        'FaceAlpha', 1, ...
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


    visualizeAMS(U0,cfg0.A,1,AMSopts,1,1,1);

    %visualizeSlider(cfg,U);
end