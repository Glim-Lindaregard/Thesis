function cfg = config()
% Define parameters 
cfg = struct();

cfg.N = 0;
cfg.beta = zeros(8,1);


cfg.u_max = 0.7*ones(8,1);      %Max thruster outputs [N]
cfg.u_min = zeros(8,1);         %Min thruster outputs [N]
a = 0.2;
cfg.pos  = [ +a,+a;
             +a,+a;
             -a,+a;
             -a,+a;
             -a,-a;
             -a,-a;
             +a,-a;
             +a,-a ];           %Thruster positions [x y] per thruster [m]
cfg.a = a;                      %Length from body center to edge [m]

cfg.N = length(cfg.pos(:,1));   %Number of thrusters

cfg.beta = [ 3*pi/2; 
                  pi; 
                   0; 
              3*pi/2; 
                pi/2; 
                   0; 
                  pi; 
                pi/2 ];         %Thruster angles from +x [radians]


% Build A matrix
A = zeros(3,cfg.N);
for i=1:cfg.N
    bx = cfg.beta(i);  c = cos(bx); s = sin(bx);
    rx = cfg.pos(i,1); ry = cfg.pos(i,2);
    A(:,i) = [ c; s; rx*s - ry*c ];
end
cfg.A = A;

end