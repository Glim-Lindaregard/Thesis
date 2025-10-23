function cfg = config()
% Define parameters for the configuration
% a: half side-length, umax: max thrust per jet [N], e.g., a=0.20; umax=0.7;
a = 0.2; 
umax = 0.7;
umin = 0;
cfg.pos  = [ +a,+a;
             +a,+a;
             -a,+a;
             -a,+a;
             -a,-a;
             -a,-a;
             +a,-a;
             +a,-a ];             % [x y] per thruster

% cfg.pos  = [ +a,+a;
%              +a,+a;
%              -a,+a;
%              -a,+a];  

m = length(cfg.pos(:,1)); %nr thrusters

%cfg.beta = [ pi/2; 0; pi/2; pi; 3*pi/2; pi; 3*pi/2; 0] + pi;   %Thruster angles from +x
cfg.beta = [ 3*pi/2; pi; 0; 3*pi/2; pi/2; 0; pi; pi/2];   %Thruster angles from +x
%cfg.beta = [ 3*pi/2; pi; 0; 3*pi/2];
cfg.u_min = umin*ones(m,1);
cfg.u_max = umax*ones(m,1);
cfg.a = a;
cfg.m = m;

% Build A
A = zeros(3,m);
for i=1:m
    bx = cfg.beta(i);  c = cos(bx); s = sin(bx);
    rx = cfg.pos(i,1); ry = cfg.pos(i,2);
    A(:,i) = [ c; s; rx*s - ry*c ];
end
cfg.A = A;  %Reverse force verctors so thr thrusters puch and not pull

end