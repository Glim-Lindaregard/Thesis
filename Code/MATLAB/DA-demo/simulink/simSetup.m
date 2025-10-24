% --- Setup sim mvariables ---


cfg = config();                  % your existing config()
AMS = buildAMS_row(cfg);         % make buildAMS_row accept cfg

% ---- basic sim parameters (adjust as you like)
simTime = 5;
simParam.StopTime  = '20';            % seconds, as string


ref = [1,1,3*pi/2,0,0,0]';

init = [0 0 pi/2 10 0 0];

% Constants
m    = 4.53; % [kg]
I_zz = 0.11; % [kgm^2]

Kp = diag([12 12 0.99]);
Kd = diag([8.7 8.7 0.53]);


% ---- inject CFG & AMS into model and set parameters
%in = Simulink.SimulationInput(simParam.model);
%in = in.setVariable('AMS', AMS);


% ---- run
simout = sim('sliderSim');


animateSliderTrajectory(simout,ref);

plotStates(simout,ref);




