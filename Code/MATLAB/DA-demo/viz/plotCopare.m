function plotCopare(t, X_no, X_yes, xRef)
% PLOTSTATESCORIOLIS
%   Compare all 6 states with and without Coriolis terms in the plant.
%
% Inputs:
%   t     : Nt x 1 or 1 x Nt time vector
%   X_no  : Nt x 6 states, plant WITHOUT Coriolis
%   X_yes : Nt x 6 states, plant WITH Coriolis
%   xRef  : 6 x 1 or 1 x 6 reference state (constant)

    % --- ensure shapes are consistent ---
    t = t(:);                    % column vector

    if size(X_no,2) ~= 6
        % assume X_no is 6 x Nt
        X_no = X_no.';           % make Nt x 6
    end
    if size(X_yes,2) ~= 6
        X_yes = X_yes.';         % make Nt x 6
    end

    xRef = xRef(:);              % 6x1

    % --- extract states ---
    x_no  = X_no(:,1);  y_no  = X_no(:,2);  th_no = X_no(:,3);
    vx_no = X_no(:,4);  vy_no = X_no(:,5);  r_no  = X_no(:,6);

    x_yes  = X_yes(:,1);  y_yes  = X_yes(:,2);  th_yes = X_yes(:,3);
    vx_yes = X_yes(:,4);  vy_yes = X_yes(:,5);  r_yes  = X_yes(:,6);

    % --- labels for plotting ---
    stateNames = {'x [m]', 'y [m]', '\theta [rad]', 'v_x [m/s]', 'v_y [m/s]', 'r [rad/s]'};

    % --- figure and subplots ---
    figure;
    
    % 1) x
    subplot(3,2,1);
    hold on; grid on;
    plot(t, x_no, 'LineWidth', 1.3);
    plot(t, x_yes, '--', 'LineWidth', 1.3);
    yline(xRef(1), ':', 'LineWidth', 1.0);
    ylabel(stateNames{1});
    title('Position and velocity comparison');
    legend('No Coriolis', 'With Coriolis', 'Reference', 'Location', 'best');

    % 2) y
    subplot(3,2,3);
    hold on; grid on;
    plot(t, y_no, 'LineWidth', 1.3);
    plot(t, y_yes, '--', 'LineWidth', 1.3);
    yline(xRef(2), ':', 'LineWidth', 1.0);
    ylabel(stateNames{2});

    % 3) theta
    subplot(3,2,5);
    hold on; grid on;
    plot(t, th_no, 'LineWidth', 1.3);
    plot(t, th_yes, '--', 'LineWidth', 1.3);
    yline(xRef(3), ':', 'LineWidth', 1.0);
    ylabel(stateNames{3});

    % 4) v_x
    subplot(3,2,2);
    hold on; grid on;
    plot(t, vx_no, 'LineWidth', 1.3);
    plot(t, vx_yes, '--', 'LineWidth', 1.3);
    yline(xRef(4), ':', 'LineWidth', 1.0);
    ylabel(stateNames{4});

    % 5) v_y
    subplot(3,2,4);
    hold on; grid on;
    plot(t, vy_no, 'LineWidth', 1.3);
    plot(t, vy_yes, '--', 'LineWidth', 1.3);
    yline(xRef(5), ':', 'LineWidth', 1.0);
    ylabel(stateNames{5});
    xlabel('t [s]');

    % 6) r
    subplot(3,2,6);
    hold on; grid on;
    plot(t, r_no, 'LineWidth', 1.3);
    plot(t, r_yes, '--', 'LineWidth', 1.3);
    yline(xRef(6), ':', 'LineWidth', 1.0);
    ylabel(stateNames{6});
    xlabel('t [s]');
end
