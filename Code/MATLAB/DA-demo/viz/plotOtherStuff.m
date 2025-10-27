function plotOtherStuff(simout)
    % -- Plotts comanded a from the PID VS Allocaded a ---


    logs = simout.logsout;

    aPD = logs.get('ad').Values;
    t_PD_T  = aPD.Time;
    aPDdata = aPD.Data;


    aAllo = logs.get('a').Values;
    t_allo_T = aAllo.Time;
    aAllodata= aAllo.Data;

    
    aPDdataX = aPDdata(1,:); aPDdataY = aPDdata(2,:); aPDdataTh = aPDdata(3,:);
    aAlloX = aAllodata(1,:); aAlloY = aAllodata(2,:); aAlloTh = aAllodata(3,:);


    t_PD = seconds(t_PD_T - t_PD_T(1));
    t_allo = seconds(t_allo_T- t_allo_T(1));


    % --- 3×1 grid of plots ---
    figure('Name','Moment tracking','Color',[0.2 0.2 0.2]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % aPDx vs aAlloX
    nexttile;
    plot(t_PD, aPDdataX, '-', 'LineWidth', 1.4, 'DisplayName','aPDx'); hold on;
    plot(t_allo, aAlloX, '-', 'LineWidth', 1.4, 'DisplayName','aAlloX'); hold on;
    grid on; ylabel('Force [N]'); legend('aPDx','aAlloX'); title('aPDx VS aAlloX');

    % aPDy vs aAlloY
    nexttile;
    plot(t_PD, aPDdataY, '-', 'LineWidth', 1.4, 'DisplayName','aPDy'); hold on;
    plot(t_allo, aAlloY, '-', 'LineWidth', 1.4, 'DisplayName','aAlloY'); hold on;
    grid on; ylabel('Force [N]'); legend('aPDy','aAlloY'); title('aPDy VS aAlloY');

    % aPDth vs aAlloTH
    nexttile;
    plot(t_PD, aPDdataTh, '-', 'LineWidth', 1.4, 'DisplayName','aPDth'); hold on;
    plot(t_allo, aAlloTh, '-', 'LineWidth', 1.4, 'DisplayName','aAlloTh'); hold on;
    grid on; ylabel('Moment [Nm]'); legend('aPDth','aAlloTh'); title('aPDth VS aAlloTh');

end