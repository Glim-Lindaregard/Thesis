% plot_Uctrl_cases_states_dark.m
% -------------------------------------------------------------------------
% Dark background + global legend.
% Computes Uctrl, removes all-healthy, plots remaining cases in one figure.
% -------------------------------------------------------------------------

clear; close all; clc;

[~, Uctrl] = d4_failure_classes_and_controllable();

% drop all-healthy
allHealthy = ones(1,8);
Uplot = Uctrl(~all(Uctrl == allHealthy, 2), :);
nCases = size(Uplot,1);
if nCases == 0
    error('No cases to plot after removing all-healthy.');
end

% grid size (near-square)
nCols = ceil(sqrt(nCases));
nRows = ceil(nCases / nCols);

cfg  = config();
a    = cfg.a;
pos  = cfg.pos(:,1:2);
beta = cfg.beta(:);
A    = cfg.A;
umax = cfg.u_max(1);

% Colors
colPassive = [0.2 0.4 0.9];   % blue
colHealthy = [0.65 0.65 0.65];% gray
colActive  = [1.0 0.55 0.1];  % orange

colFig   = [0.08 0.08 0.10];  % dark background
colAx    = [0.08 0.08 0.10];  % same for axes
colEdge  = [0.70 0.70 0.75];  % light edges
colText  = [0.92 0.92 0.95];  % near-white text

% --- Figure
fig = figure('Name','Uctrl cases (state: passive/healthy/active)', ...
    'Color', colFig);
sgtitle(sprintf('Controllable D4-unique cases (excluding all-healthy): %d', nCases), ...
    'Color', colText);

for k = 1:nCases
    s = Uplot(k,:); % 0/1/2

    ax = subplot(nRows, nCols, k);
    set(ax, 'Color', colAx);
    cla(ax); hold(ax,'on'); axis(ax,'equal');

    xlim(ax, [-1.2*a 1.2*a]); ylim(ax, [-1.2*a 1.2*a]);

    % table boundary (light on dark)
    rectangle(ax, 'Position',[-a -a 2*a 2*a], 'EdgeColor', colEdge);

    % thruster markers by state
    for i = 1:8
        if s(i) == 0
            plot(ax, pos(i,1), pos(i,2), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', colPassive, 'MarkerEdgeColor', colEdge);
        elseif s(i) == 1
            plot(ax, pos(i,1), pos(i,2), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', colHealthy, 'MarkerEdgeColor', colEdge);
        elseif s(i) == 2
            plot(ax, pos(i,1), pos(i,2), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', colActive, 'MarkerEdgeColor', colEdge);
        else
            error('Unexpected state value (must be 0/1/2).');
        end
    end

    % stuck-on direction arrows (only state==2)
    Lmax = 0.45*a;
    for i = 1:8
        if s(i) == 2
            d = [cos(beta(i)) sin(beta(i))];
            quiver(ax, pos(i,1), pos(i,2), Lmax*d(1), Lmax*d(2), 0, ...
                'LineWidth', 1.2, 'MaxHeadSize', 0.8, 'Color', colActive);
        end
    end

    % visualize resultant wrench from stuck-on thrusters only
    ud = zeros(8,1);
    ud(s == 2) = umax;

    y = A * ud;
    f = y(1:2);
    tau = y(3);

    % resultant force (red)
    sf = (0.6*a) / max(norm(f), eps);
    quiver(ax, 0,0, sf*f(1), sf*f(2), 0, 'LineWidth',2, ...
        'Color',[1 0 0], 'MaxHeadSize',1.0);

    % torque arc (magenta)
    R = 0.7*a;
    tau_scale = a * sum(abs(ud));
    ang = min(2*pi, abs(tau)/max(tau_scale,eps)*2*pi);
    if ang > 0
        sgn = sign(tau); if sgn==0, sgn=1; end
        t = linspace(0, sgn*ang, 100);
        xx = R*cos(t); yy = R*sin(t);
        plot(ax, xx, yy, '-', 'LineWidth',1.6, 'Color',[1 0 1]);

        te = sgn*ang;
        tv = R*[-sin(te); cos(te)]; tv = tv / max(norm(tv),eps);
        quiver(ax, R*cos(te), R*sin(te), 0.1*a*tv(1), 0.1*a*tv(2), 0, ...
            'Color',[1 0 1], 'MaxHeadSize',1.2);
    end

    % subplot label
    title(ax, sprintf('%d', k), 'Color', colText, 'FontSize', 9);
    axis(ax,'off');
    hold(ax,'off');
end

% ----------------- Global Legend (single) -----------------
% Create a tiny invisible axes and draw "dummy" objects to feed legend().
axL = axes('Position',[0.01 0.01 0.01 0.01], 'Visible','off', 'Color','none');
hold(axL,'on');

hP = plot(axL, nan,nan,'o','MarkerSize',7,'MarkerFaceColor',colPassive,'MarkerEdgeColor',colEdge);
hH = plot(axL, nan,nan,'o','MarkerSize',7,'MarkerFaceColor',colHealthy,'MarkerEdgeColor',colEdge);
hA = plot(axL, nan,nan,'o','MarkerSize',7,'MarkerFaceColor',colActive, 'MarkerEdgeColor',colEdge);

hActArrow = quiver(axL, nan,nan, 1,0, 0, 'LineWidth',1.2, 'Color',colActive, 'MaxHeadSize',0.8);
hF       = quiver(axL, nan,nan, 1,0, 0, 'LineWidth',2.0, 'Color',[1 0 0], 'MaxHeadSize',1.0);
hTau     = plot(axL, nan,nan, '-', 'LineWidth',1.6, 'Color',[1 0 1]);

leg = legend(axL, [hP hH hA hActArrow hF hTau], ...
    {'Passive failed (0)', 'Healthy/Neutral (1)', 'Active stuck-on (2)', ...
     'Stuck-on direction', 'Resultant force (from stuck-on)', 'Resultant torque (from stuck-on)'}, ...
    'TextColor', colText, 'Color', 'none', 'EdgeColor', colEdge);

% Put legend at bottom; tweak if you want it elsewhere
%set(leg, 'Location','southoutside', 'Orientation','horizontal', 'Box','on');

% Make legend font a bit smaller to fit
%leg.FontSize = 9;

hold(axL,'off');
