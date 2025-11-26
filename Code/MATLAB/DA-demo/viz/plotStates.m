function plotStates_arrays(t, X, ref)
% PLOTSTATES_ARRAYS  Plot state tracking using array-based logs.
%   t   : 1xN or Nx1 time vector [s]
%   X   : 6xN state trajectory [x; y; theta; vx; vy; w]
%   ref : 6x1 reference state [x_ref; y_ref; theta_ref; vx_ref; vy_ref; w_ref]

    % Ensure t is a column and starts at 0
    t = t(:);
    t = t - t(1);

    % X is 6xN; transpose to N×6 for column-wise access
    X = X.';
    x  = X(:,1);
    y  = X(:,2);
    th = X(:,3);
    vx = X(:,4);
    vy = X(:,5);
    w  = X(:,6);

    % With your new model, vx, vy are already world-frame velocities
    xdot = vx;
    ydot = vy;

    % Constant reference values (broadcast ref over time)
    xref    = ref(1) * ones(size(x));
    yref    = ref(2) * ones(size(y));
    thref   = ref(3) * ones(size(th));
    xdotref = ref(4) * ones(size(xdot));
    ydotref = ref(5) * ones(size(ydot));
    wref    = ref(6) * ones(size(w));

    % --- figure & layout ---
    fig = figure('Name','State tracking','Color',[0.2 0.2 0.2]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % x vs xref
    nexttile;
    plot(t, x, '-', 'LineWidth', 1.4, 'DisplayName','$x$'); hold on;
    plot(t, xref, '--', 'LineWidth', 1.2, 'DisplayName','$x_{\mathrm{ref}}$');
    grid on;
    ylabel('$x$ [m]','Interpreter','latex');
    title('$x$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');

    % xdot vs xdotref
    nexttile;
    plot(t, xdot, '-', 'LineWidth', 1.4, 'DisplayName','$\dot{x}$'); hold on;
    plot(t, xdotref, '--', 'LineWidth', 1.2, 'DisplayName','$\dot{x}_{\mathrm{ref}}$');
    grid on;
    ylabel('$\dot{x}$ [m/s]','Interpreter','latex');
    title('$\dot{x}$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');

    % y vs yref
    nexttile;
    plot(t, y, '-', 'LineWidth', 1.4, 'DisplayName','$y$'); hold on;
    plot(t, yref, '--', 'LineWidth', 1.2, 'DisplayName','$y_{\mathrm{ref}}$');
    grid on;
    ylabel('$y$ [m]','Interpreter','latex');
    title('$y$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');

    % ydot vs ydotref
    nexttile;
    plot(t, ydot, '-', 'LineWidth', 1.4, 'DisplayName','$\dot{y}$'); hold on;
    plot(t, ydotref, '--', 'LineWidth', 1.2, 'DisplayName','$\dot{y}_{\mathrm{ref}}$');
    grid on;
    ylabel('$\dot{y}$ [m/s]','Interpreter','latex');
    title('$\dot{y}$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');

    % theta vs thetaref
    nexttile;
    plot(t, th, '-', 'LineWidth', 1.4, 'DisplayName','$\theta$'); hold on;
    plot(t, thref, '--', 'LineWidth', 1.2, 'DisplayName','$\theta_{\mathrm{ref}}$');
    grid on;
    ylabel('$\theta$ [rad]','Interpreter','latex');
    xlabel('Time [s]','Interpreter','latex');
    title('$\theta$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');

    % omega vs omegaref
    nexttile;
    plot(t, w, '-', 'LineWidth', 1.4, 'DisplayName','$\omega$'); hold on;
    plot(t, wref, '--', 'LineWidth', 1.2, 'DisplayName','$\omega_{\mathrm{ref}}$');
    grid on;
    ylabel('$\omega$ [rad/s]','Interpreter','latex');
    xlabel('Time [s]','Interpreter','latex');
    title('$\omega$ tracking','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
end
