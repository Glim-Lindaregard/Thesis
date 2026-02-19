%% Lissajous circle (single plot)
clear; clc;

% Parameters (same meaning as ROS)
A     = 0.5;    % lissa_A
B     = 0.5;    % lissa_B
a     = 1;      % lissa_a
b     = 0;      % lissa_b
delta = 0;    % lissa_delta
omega = 0.2;    % lissa_omega

T = 2*pi/omega;          % one full cycle
t = linspace(0, T, 2000);

% Lissajous equations (identical to Python)
theta = omega * t;
x = A * sin(a * theta + delta);
y = B * sin(b * theta);

% Plot
figure;
plot(x, y, 'LineWidth', 2);
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('Lissajous trajectory (circle)');
