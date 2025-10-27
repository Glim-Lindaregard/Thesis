function plotStates(simout, ref)

    % --- Extract states ---
    logs = simout.logsout;
    State = logs.get('State').Values;
    t_state  = State.Time;
    X        = State.Data;


    x = X(:,1); y = X(:,2); th = X(:,3);
    vx = X(:,4);  vy = X(:,5);  w = X(:,6);

    
    ct = cos(th); st = sin(th);
    
    % rotate body-frame velocities to world frame
    xdot =  ct.*vx - st.*vy;
    ydot =  st.*vx + ct.*vy;


    t = seconds(t_state - t_state(1));   % time vector starting at 0 s

    % --- Constant reference values ---
    xref  = ref(1) * ones(size(x));
    yref  = ref(2) * ones(size(y));
    thref = ref(3) * ones(size(th));
    xdotref = ref(4) * ones(size(xdot));
    ydotref = ref(5) * ones(size(ydot));
    wref  = ref(6) * ones(size(w));

    % --- 3×2 grid of plots ---
    figure('Name','State tracking','Color',[0.2 0.2 0.2]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % x vs xref
    nexttile;
    plot(t, x, '-', 'LineWidth', 1.4, 'DisplayName','x'); hold on;
    plot(t, xref, '--', 'LineWidth', 1.2, 'DisplayName','x_{ref}');
    grid on; ylabel('x [m]'); legend('Location','best'); title('x tracking');

    % xdot vs xref
    nexttile;
    plot(t, xdot, '-', 'LineWidth', 1.4, 'DisplayName','x velocity'); hold on;
    plot(t, xdotref, '--', 'LineWidth', 1.2, 'DisplayName','xdot_{ref}');
    grid on; ylabel('xdot [m/s]'); legend('x velocity','xdot_{ref}'); title('xdot tracking');

    % y vs yref
    nexttile;
    plot(t, y, '-', 'LineWidth', 1.4, 'DisplayName','y'); hold on;
    plot(t, yref, '--', 'LineWidth', 1.2, 'DisplayName','y_{ref}');
    grid on; ylabel('y [m]'); legend('Location','best'); title('y tracking');

    % ydot vs yref
    nexttile;
    plot(t, ydot, '-', 'LineWidth', 1.4, 'DisplayName','y velocity'); hold on;
    plot(t, ydotref, '--', 'LineWidth', 1.2, 'DisplayName','ydot_{ref}');
    grid on; ylabel('ydot [m/s]'); legend('y velocity','ydot_{ref}'); title('ydot tracking');

    % theta vs thetaref
    nexttile;
    plot(t, th, '-', 'LineWidth', 1.4, 'DisplayName','\theta'); hold on;
    plot(t, thref, '--', 'LineWidth', 1.2, 'DisplayName','\theta_{ref}');
    grid on; ylabel('\theta [rad]'); xlabel('Time [s]');
    legend('Location','best'); title('\theta tracking');

    % omega vs thetaref
    nexttile;
    plot(t, w, '-', 'LineWidth', 1.4, 'DisplayName','\omega'); hold on;
    plot(t, wref, '--', 'LineWidth', 1.2, 'DisplayName','\omega_{ref}');
    grid on; ylabel('\omega [rad/s]'); xlabel('Time [s]');
    legend('\omega','\omega_{ref}'); title('\omega tracking');

end
