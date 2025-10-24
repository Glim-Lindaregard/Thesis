function plotStates(simout, ref)
% PLOTSLIDERTRAJECTORY  Plot x, y, theta vs constant reference values
%
% Usage:
%   plotSliderTrajectory(simout, [x_ref; y_ref; theta_ref])
%
% Extracts pose data from simout.logsout ('state') using the same
% logic as animateSliderTrajectory, then plots x/xref, y/yref,
% and theta/theta_ref in a 3×1 grid.

    % --- Extract pose (same as your animate function) ---
    logs = simout.logsout;
    state_ts = logs.get('state').Values;
    t_state  = state_ts.Time;
    X        = state_ts.Data;    % could be N×3 or 3×N depending on block

    if size(X,2)==3
        x = X(:,1); y = X(:,2); th = X(:,3);
    elseif size(X,1)==3
        x = X(1,:).'; y = X(2,:).'; th = X(3,:).';
    else
        error('Logged ''state'' must have width 3. Got size %s.', mat2str(size(X)));
    end

    t = seconds(t_state - t_state(1));   % time vector starting at 0 s

    % --- Constant reference values ---
    xref  = ref(1) * ones(size(x));
    yref  = ref(2) * ones(size(y));
    thref = ref(3) * ones(size(th));

    % --- 3×1 grid of plots ---
    figure('Name','States tracking','Color',[0.2 0.2 0.2]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % x vs xref
    nexttile;
    plot(t, x, '-', 'LineWidth', 1.4, 'DisplayName','x'); hold on;
    plot(t, xref, '--', 'LineWidth', 1.2, 'DisplayName','x_{ref}');
    grid on; ylabel('x [m]'); legend('Location','best'); title('x tracking');

    % y vs yref
    nexttile;
    plot(t, y, '-', 'LineWidth', 1.4, 'DisplayName','y'); hold on;
    plot(t, yref, '--', 'LineWidth', 1.2, 'DisplayName','y_{ref}');
    grid on; ylabel('y [m]'); legend('Location','best'); title('y tracking');

    % theta vs thetaref
    nexttile;
    plot(t, th, '-', 'LineWidth', 1.4, 'DisplayName','\theta'); hold on;
    plot(t, thref, '--', 'LineWidth', 1.2, 'DisplayName','\theta_{ref}');
    grid on; ylabel('\theta [rad]'); xlabel('Time [s]');
    legend('Location','best'); title('\theta tracking');

end
