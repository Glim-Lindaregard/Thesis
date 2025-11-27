function cfg = config()
cfg = struct();

% Simulation parameters 
cfg.Ts       = 0.2;
cfg.simTime = 20;
cfg.x0   = [0; 0; pi/2; 0; 0; 0];
cfg.xRef = [2; 0; pi/2; 0; 0; 0];

%MPC parameters
cfg.FxMax   = 0.8*(2*0.7);        % tune
cfg.FyMax   = 0.8*(2*0.7);
cfg.TauMax  = 0.8*(4*0.14);

cfg.N = 15;

cfg.Q = diag([2 2 2  10 10 10]);   % << try this

cfg.R = diag(1*[2 2 0.8]);          % big damping increase


%Physical
cfg.m    = 4.436; % [kg]
cfg.Izz = 1.092; % [kgm^2]

cfg.u_max = 0.7*ones(8,1);      %Max thruster outputs [N]
cfg.u_min = zeros(8,1);         %Min thruster outputs [N]

a = 0.195; b = 0.140;

cfg.pos = [ +a,+b;
            +b,+a;
            -b,+a;   % fixed
            -a,+b;
            -a,-b;
            -b,-a;
            +b,-a;
            +a,-b ];

cfg.a = a;                      %Length from body center to edge [m]

cfg.N_thrusters = length(cfg.pos(:,1));   %Number of thrusters

cfg.beta = [0, pi/2, pi/2, pi, pi, 3*pi/2, 3*pi/2, 0]' + pi;



% Build A matrix
A = zeros(3,cfg.N_thrusters);
for i=1:cfg.N_thrusters
    bx = cfg.beta(i);  c = cos(bx); s = sin(bx);
    rx = cfg.pos(i,1); ry = cfg.pos(i,2);
    A(:,i) = [ c; s; rx*s - ry*c ];
end
cfg.A = A;

end