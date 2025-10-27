close all;
% --- Setup sim variables ---
cfg = config();
AMS = buildAMS_row(cfg);

simTime = 25;

ref = [2,1,-pi/2,0,0,0]';

init = [0 0 pi/2 0 0 0];

% Constants
m    = 4.436; % [kg]
I_zz = 1.092; % [kgm^2]


% PID parameters (TEMP)
Kp = diag([1 1 0.48]);
Kd = diag([2 2 1.2]);
ki = diag([0.06 0.06 0.064]);


% --- run ---
simout = sim('sliderSim');


%Visualizion
opts.fps = 250;
opts.saveVideo = false;
opts.videoName = 'sliderAnimation.mp4';
opts.simSpeed = 60;
animateTrajectory(simout,ref,opts);

%plotStates(simout,ref);

%plotOtherStuff(simout);




