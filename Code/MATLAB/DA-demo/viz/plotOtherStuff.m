function plotOtherStuff(t, adHist, aRealHist)
% PLOTOTHERSTUFF  Compare commanded wrench vs allocated wrench
%
% Inputs:
%   t         : 1xNt time vector (same as simulateSlider output)
%   adHist    : 3x(Nt-1) commanded wrench from MPC
%   aRealHist : 3x(Nt-1) realized wrench after allocation

    % Time stamps for a_d and a_real (they live on k=1..Nt-1)
    t_cmd = t(1:end-1);

    % Split components
    aPDdataX   = adHist(1,:);   % commanded Fx
    aPDdataY   = adHist(2,:);   % commanded Fy
    aPDdataTh  = adHist(3,:);   % commanded tau

    aAlloX     = aRealHist(1,:);  % realized Fx
    aAlloY     = aRealHist(2,:);  % realized Fy
    aAlloTh    = aRealHist(3,:);  % realized tau

    % --- 3×1 grid of plots ---
    figure('Name','Wrench tracking','Color',[0.2 0.2 0.2]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % Fx
    nexttile;
    plot(t_cmd, aPDdataX, '-',  'LineWidth', 1.4); hold on;
    plot(t_cmd, aAlloX,  '--', 'LineWidth', 1.4);
    grid on;
    ylabel('Force [N]');
    legend('a_d,x','a_{real,x}','Location','best');
    title('Commanded vs realized F_x','Interpreter','tex');

    % Fy
    nexttile;
    plot(t_cmd, aPDdataY, '-',  'LineWidth', 1.4); hold on;
    plot(t_cmd, aAlloY,  '--', 'LineWidth', 1.4);
    grid on;
    ylabel('Force [N]');
    legend('a_d,y','a_{real,y}','Location','best');
    title('Commanded vs realized F_y','Interpreter','tex');

    % Tau
    nexttile;
    plot(t_cmd, aPDdataTh, '-',  'LineWidth', 1.4); hold on;
    plot(t_cmd, aAlloTh,   '--', 'LineWidth', 1.4);
    grid on;
    ylabel('Moment [Nm]');
    xlabel('Time [s]');
    legend('a_d,\theta','a_{real,\theta}','Location','best');
    title('Commanded vs realized \tau','Interpreter','tex');
end
