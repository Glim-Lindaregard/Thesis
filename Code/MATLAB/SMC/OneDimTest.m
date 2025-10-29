close all

k1 = 3;  k2 = 4;

p = 2;
c = 1.5;
sigma = @(x) x(2)*c + x(1);

f = @(x,t) sin(2*t);

% Control law
epsi = 0.1;
mysgn = @(s) s / (abs(s) + epsi);

u = @(x) -c*x(2) -p*mysgn(sigma(x));

% Dynamics
f_sys = @(t,x) [ x(2);
                 u(x) + f(x,t) ];

% Initial state
x0 = [1; -2]; 


[t,X] = ode45(f_sys,[0 10],x0);

figure; 

plot(t, X(:,1), '-','Color', 'b');
hold on
plot(t,X(:,2),'--','Color','b');
sigma_vals = arrayfun(@(i) sigma(X(i,:)'), 1:length(t));
plot(t, sigma_vals);
grid on;

figure;
plot(X(:,1), X(:,2));

figure;
u_vals = arrayfun(@(i) u(X(i,:)'), 1:length(t));
plot(t,u_vals);

figure; 
fplot(mysgn,[-5 5]);

