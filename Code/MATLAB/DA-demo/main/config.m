function cfg = config()
cfg = struct();

% Simulation parameters 
cfg.useCoriolis = 1;
cfg.Ts       = 0.2;
cfg.simTime = 15;

%MPC parameters
cfg.FxMax   = 1.5*(2*0.7);        % tune
cfg.FyMax   = 1.5*(2*0.7);
cfg.TauMax  = 1.5*(4*0.14);

cfg.N = 15;

cfg.Q = diag([10 10 10  2 2 2]);   % << try this

cfg.R = diag(0.1*[2 2 0.8]);          % big damping increase


%Physical
cfg.m    = 4.436; % [kg]
cfg.Izz = 1.092; % [kgm^2]
a = 0.195; b = 0.140;
cfg.a = a;                      %Length from body center to edge [m]

wallBuffer = 2*cfg.a;
cfg.xMin = wallBuffer;
cfg.xMax = 4 - wallBuffer;
cfg.yMin = wallBuffer;
cfg.yMax = 4 - wallBuffer;
cfg.wallBuffer = wallBuffer;

cfg.u_max = 0.7*ones(8,1);      %Max thruster outputs [N]
cfg.u_min = zeros(8,1);         %Min thruster outputs [N]


%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%
cfg.x0   = [1; 1;0 ;0;0;0];
cfg.xRef = [2; 2; pi; 0; 0; 0]; 

r     = 1.0;
omega = 0.2;  % rad/s

% cfg.refFun = @(t) [ ...
%     1.5 + r*cos(omega*t);         % x(t)
%     1.5 + r*sin(omega*t);         % y(t)
%     atan2(sin(omega*t),cos(omega*t));  % theta ~ pointing along the path
%     0;
%     0;
%     0];

% cfg.refFun = @(t) [1;1;0;
%     0;
%     0;
%     0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.pos = [ +a,+b;
            +b,+a;
            -b,+a;   % fixed
            -a,+b;
            -a,-b;
            -b,-a;
            +b,-a;
            +a,-b ];


cfg.N_thrusters = length(cfg.pos(:,1));   %Number of thrusters

cfg.beta = [0, pi/2, pi/2, pi, pi, 3*pi/2, 3*pi/2, 0]' + pi;



% Build A matrix
A = zeros(3,cfg.N_thrusters);
for i=1:cfg.N_thrusters
    bx = cfg.beta(i);  c = cos(bx); s = sin(bx);
    rx = cfg.pos(i,1); ry = cfg.pos(i,2);
    A(:,i) = [ c; s; rx*s - ry*c ];
end

A =  [-1,     0,     1,     0,     1,     0,    -1,     0;
       0,     1,     0,     1,     0,    -1,     0,    -1;
   -0.14,  0.14,  0.14, -0.14, -0.14,  0.14,  0.14, -0.14];

cfg.A = A;

end