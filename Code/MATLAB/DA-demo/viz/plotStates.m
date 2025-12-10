function plotStates(t, X, ref)
% PLOTSTATES_ARRAYS  Plot state tracking using array-based logs.
%   t   : 1xN or Nx1 time vector [s]
%   X   : 6xN state trajectory [x; y; theta; vx; vy; w]
%   ref : one of
%           - 6x1 constant reference
%           - 6xN time-varying reference
%           - function handle: ref(tk) -> 6x1 at time tk

    % Ensure t is a column and starts at 0
    t = t(:);
    t = t - t(1);
    N = numel(t);

    % X is 6xN; transpose to N×6 for column-wise access
    if size(X,1) ~= 6 && size(X,2) == 6
        X = X.';   % allow Nx6 input
    elseif size(X,1) == 6
        X = X.';   % 6xN -> Nx6
    else
        error('X must be 6xN or Nx6.');
    end

    if size(X,1) ~= N
        error('Time vector length and state trajectory length must match.');
    end

    x  = X(:,1);
    y  = X(:,2);
    th = X(:,3);
    vx = X(:,4);
    vy = X(:,5);
    w  = X(:,6);

    % With your new model, vx, vy are already world-frame velocities
    xdot = vx;
    ydot = vy;

    % ---- build reference trajectories in array form: N×6 ----
    if isa(ref, 'function_handle')
        % time-varying via function handle
        Ref = zeros(N,6);
        for k = 1:N
            rk = ref(t(k));          % expect 6x1
            rk = rk(:);
            if numel(rk) ~= 6
                error('ref(t) must return a 6x1 vector.');
            end
            Ref(k,:) = rk.';
        end

    elseif isvector(ref) && numel(ref) == 6
        % constant reference (old behaviour)
        r = ref(:).';
        Ref = repmat(r, N, 1);       % N×6

    else
        % array reference: expect 6xN or Nx6
        if size(ref,1) == 6
            Ref = ref.';             % 6xN -> N×6
        elseif size(ref,2) == 6
            Ref = ref;               % Nx6 already
        else
            error('ref must be 6x1, 6xN, Nx6, or a function handle.');
        end

        if size(Ref,1) ~= N
            error('ref length must match time vector length.');
        end
    end

    % Split reference into components
    xref    = Ref(:,1);
    yref    = Ref(:,2);
    thref   = Ref(:,3);
    xdotref = Ref(:,4);
    ydotref = Ref(:,5);
    wref    = Ref(:,6);

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
