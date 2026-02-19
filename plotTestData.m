% pwm_sweep_plot.m
% Read sweep results (frequency, resolution, clicks, x_err, t_err)
% NOTE: CSV columns still named x_err and t_err, but they now represent:
%   x_err -> RMSE_pos
%   t_err -> ITAE_pos
%
% and generate:
%  1) clicks vs resolution (one figure, multiple frequency traces)
%  2) RMSE_pos vs resolution (one figure, multiple frequency traces)
%  3) ITAE_pos vs resolution (one figure, multiple frequency traces)
%  4) clicks vs frequency (one figure, multiple resolution traces)
%  5) RMSE_pos vs frequency (one figure, multiple resolution traces)
%  6) ITAE_pos vs frequency (one figure, multiple resolution traces)
%  7) heatmaps with SQUARE cells (categorical spacing)
%     - white figure surround, white plot background
%     - ALL text/labels/axes black
%     - colorbar tick numbers black + evenly spaced ticks

clear; clc; close all;

% ------------------ Config ------------------
csvFile = "sweep_results_V2.csv";   % <- change if needed

% ------------------ Load --------------------
T = readtable(csvFile);

requiredCols = ["frequency","resolution","clicks","x_err","t_err"];  % CSV still uses these names
missing = setdiff(requiredCols, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("CSV is missing columns: %s", strjoin(missing, ", "));
end

% Ensure numeric (Google Sheets sometimes exports as text)
T.frequency   = double(T.frequency);
T.resolution  = double(T.resolution);
T.clicks      = double(T.clicks./50);
T.x_err       = double(T.x_err);   % now RMSE_pos
T.t_err       = double(T.t_err);   % now ITAE_pos

% Sort for nicer plots
T = sortrows(T, ["frequency","resolution"]);

freqs   = unique(T.frequency);
resVals = unique(T.resolution);

%------------------- cost function ----------------------
% Weights (tune these)
wRMSE   = 0.35;
wITAE   = 0.35;
wClicks = 0.30;

best = findBestFreqRes(T, wRMSE, wITAE, wClicks, 95, 2.0);

disp("=== Best (frequency,resolution) ===");
disp(best.table(1:10, :));   % top 10 candidates
fprintf("BEST: f=%g, res=%g | cost=%.3f | clicks=%.1f | RMSE=%.4f | ITAE=%.4f\n", ...
    best.frequency, best.resolution, best.cost, best.clicks, best.RMSE_pos, best.ITAE_pos);

% ------------------ Multi-trace plots ------------------
% metric vs resolution (traces = frequency)
plotMetricVsResolution(T, freqs, "clicks", "Clicks/s", ...
    "Clicks/s vs Resolution (grouped by Frequency)");

plotMetricVsResolution(T, freqs, "x_err", "RMSE_{pos}", ...
    "RMSE_{pos} vs Resolution (grouped by Frequency)");

plotMetricVsResolution(T, freqs, "t_err", "ITAE_{pos}", ...
    "ITAE_{pos} vs Resolution (grouped by Frequency)");

% metric vs frequency (traces = resolution)
plotMetricVsFrequency(T, resVals, "clicks", "Clicks/s", ...
    "Clicks/s vs Frequency (grouped by Resolution)");

plotMetricVsFrequency(T, resVals, "x_err", "RMSE_{pos}", ...
    "RMSE_{pos} vs Frequency (grouped by Resolution)");

plotMetricVsFrequency(T, resVals, "t_err", "ITAE_{pos}", ...
    "ITAE_{pos} vs Frequency (grouped by Resolution)");

% ------------------ Build heatmap matrices (freq x res) ------------------
Z_clicks = nan(length(freqs), length(resVals));
Z_rmse   = nan(length(freqs), length(resVals));  % from x_err
Z_itae   = nan(length(freqs), length(resVals));  % from t_err

for i = 1:height(T)
    fi = find(freqs == T.frequency(i), 1);
    ri = find(resVals == T.resolution(i), 1);

    Z_clicks(fi, ri) = T.clicks(i);
    Z_rmse(fi, ri)   = T.x_err(i);   % RMSE_pos
    Z_itae(fi, ri)   = T.t_err(i);   % ITAE_pos
end

% ------------------ Heatmaps (SQUARE CELLS) ------------------
plotHeatmapCategorical(Z_clicks, resVals, freqs, "Resolution", "Frequency", "Clicks/s heatmap", false);
plotHeatmapCategorical(Z_rmse,   resVals, freqs, "Resolution", "Frequency", "RMSE_{pos} heatmap", true);
plotHeatmapCategorical(Z_itae,   resVals, freqs, "Resolution", "Frequency", "ITAE_{pos} heatmap",true);

%-------------------3D Heatmaps -------------------------------
plot3DSweep(resVals, freqs, Z_rmse, Z_clicks, ...
    "RMSE_{pos}", "Clicks", ...
    "3D Sweep: RMSE_{pos} (height) + Clicks/s (color)");

plot3DSweep(resVals, freqs, Z_itae, Z_clicks, ...
    "ITAE_{pos}", "Clicks", ...
    "3D Sweep: ITAE_{pos} (height) + Clicks/s (color)");

% ================== Local functions ==================

function plotMetricVsResolution(T, freqs, metricName, yLabelText, figTitle)
    figure('Color','w');
    ax = axes; %#ok<LAXES>
    set(ax, 'Color','w');
    hold(ax, 'on'); grid(ax, 'on');
    set(ax,'FontSize', 12)

    for k = 1:numel(freqs)
        f = freqs(k);
        Tf = T(T.frequency == f, :);

        % If duplicate (frequency,resolution) rows exist, average them
        [G, res] = findgroups(Tf.resolution);
        y = splitapply(@mean, Tf.(metricName), G);

        plot(ax, res, y, '-o', 'DisplayName', sprintf('f = %g', f), 'LineWidth',1.4);
    end

    xlabel(ax, 'Resolution', 'Color','k', 'FontSize',16);
    ylabel(ax, yLabelText,   'Color','k', 'FontSize',16);
    title(ax, figTitle,      'Color','k', 'FontSize',26);

    leg = legend(ax, 'Location','best');
    set(leg, 'TextColor','k', 'Color','w');

    % Force ALL axis text black
    set(ax, 'XColor','k', 'YColor','k', 'ZColor','k');
    ax.XAxis.Color = 'k';
    ax.YAxis.Color = 'k';
end

function plotMetricVsFrequency(T, resVals, metricName, yLabelText, figTitle)
    figure('Color','w');
    ax = axes; %#ok<LAXES>
    set(ax, 'Color','w');
    hold(ax, 'on'); grid(ax, 'on');
    set(ax,'FontSize', 12)

    for k = 1:numel(resVals)
        r = resVals(k);
        Tr = T(T.resolution == r, :);

        % If duplicate (frequency,resolution) rows exist, average them
        [G, f] = findgroups(Tr.frequency);
        y = splitapply(@mean, Tr.(metricName), G);

        plot(ax, f, y, '-o', 'DisplayName', sprintf('res = %g', r), 'LineWidth',1.4);
    end

    xlabel(ax, 'Frequency', 'Color','k', 'FontSize',16);
    ylabel(ax, yLabelText,  'Color','k','FontSize',16);
    title(ax, figTitle,     'Color','k','FontSize',26);

    leg = legend(ax, 'Location','best');
    set(leg, 'TextColor','k', 'Color','w','LineWidth',1.2);

    set(ax, 'XColor','k', 'YColor','k', 'ZColor','k');
    ax.XAxis.Color = 'k';
    ax.YAxis.Color = 'k';
end

function plotHeatmapCategorical(Z, resVals, freqVals, xLabelText, yLabelText, figTitle, Log)
    figure('Color','w');
    ax = axes; %#ok<LAXES>
    set(ax, 'Color','w');
    set(ax,'FontSize', 12)



    %if Log
    %     Z = log10(1+log10(1+Z));
    % end
    if Log
     Z(1,1) = NaN;
    end


    imagesc(ax, Z);

    axis(ax, 'equal');
    axis(ax, 'tight');
    set(ax, 'YDir', 'normal'); % frequency increases upward

    colormap(ax, nebula);

    xlabel(ax, xLabelText, 'Color','k', 'FontSize',16);
    ylabel(ax, yLabelText, 'Color','k','FontSize',16);
    title(ax, figTitle,    'Color','k','FontSize',26);

    % Categorical spacing: ticks are indices, labels are actual values
    xticks(ax, 1:numel(resVals));
    xticklabels(ax, string(resVals));

    yticks(ax, 1:numel(freqVals));
    yticklabels(ax, string(freqVals));

    % Make axes/ticks black (on white plot area)
    set(ax, 'XColor','k', 'YColor','k', 'ZColor','k');
    ax.XAxis.Color = 'k';
    ax.YAxis.Color = 'k';

    % Colorbar: black numbers, evenly spaced ticks
    cb = colorbar(ax);
    cb.Color = 'k';

    zMin = min(Z(:), [], 'omitnan');
    zMax = max(Z(:), [], 'omitnan');

    if isfinite(zMin) && isfinite(zMax) && zMax > zMin
        nTicks = 5;
        ticks = linspace(zMin, zMax, nTicks);
        cb.Ticks = ticks;
        cb.TickLabels = compose('%.2f', ticks);
    end
end
function plot3DSweep(resVals, freqVals, Z_height, Z_color, ...
                     zLabelText, colorLabelText, figTitle)

    figure('Color','w');
    ax = axes; hold(ax,'on'); grid(ax,'on');
    set(ax,'Color','w');
    ax.FontSize = 14;

    % Build grid
    [RES, FREQ] = meshgrid(resVals, freqVals);

    % Surface: height = performance, color = clicks
    surf(ax, RES, FREQ, Z_height, Z_color, ...
        'EdgeColor','none');

    xlabel(ax,'Resolution','Color','k');
    ylabel(ax,'Frequency','Color','k');
    zlabel(ax,zLabelText,'Color','k');
    title(ax,figTitle,'Color','k');

    set(ax,'XColor','k','YColor','k','ZColor','k');
    view(45,30);

    colormap(ax, turbo);
    cb = colorbar(ax);
    cb.Color = 'k';
    cb.Label.String = colorLabelText;
    cb.Label.Color = 'k';

    shading interp;
end
function best = findBestFreqRes(T, wRMSE, wITAE, wClicks, prc, cap)
%FINDBESTFREQRES  Pick best (frequency,resolution) using weighted normalized cost.
%   T must contain: frequency, resolution, clicks, x_err (RMSE_pos), t_err (ITAE_pos)
%   wRMSE,wITAE,wClicks: weights (e.g., 0.5,0.3,0.2)
%   prc: percentile for robust scaling (e.g., 95)
%   cap: cap normalized values (e.g., 2.0). Use Inf to disable.

    if nargin < 5 || isempty(prc), prc = 95; end
    if nargin < 6 || isempty(cap), cap = 2.0; end

    % Aggregate duplicates by mean (same freq,res)
    [G, F, R] = findgroups(T.frequency, T.resolution);
    clicks = splitapply(@mean, T.clicks, G);
    rmse   = splitapply(@mean, T.x_err,  G);   % RMSE_pos
    itae   = splitapply(@mean, T.t_err,  G);   % ITAE_pos

    % Robust scales
    sC = prctile(clicks, prc);
    sR = prctile(rmse,   prc);
    sI = prctile(itae,   prc);

    % Avoid divide-by-zero
    sC = max(sC, eps); sR = max(sR, eps); sI = max(sI, eps);

    % Normalized (robust) + optional cap
    Cn = clicks / sC;  Rn = rmse / sR;  In = itae / sI;
    if isfinite(cap)
        Cn = min(Cn, cap); Rn = min(Rn, cap); In = min(In, cap);
    end

    % Cost
    J = wRMSE*Rn + wITAE*In + wClicks*Cn;

    % Best
    [Jbest, idx] = min(J);

    best.frequency  = F(idx);
    best.resolution = R(idx);
    best.cost       = Jbest;

    best.clicks     = clicks(idx);
    best.RMSE_pos   = rmse(idx);
    best.ITAE_pos   = itae(idx);

    % also return table for inspection/sorting
    best.table = table(F, R, clicks, rmse, itae, Cn, Rn, In, J, ...
        'VariableNames', {'frequency','resolution','clicks','RMSE_pos','ITAE_pos', ...
                          'clicks_n','RMSE_n','ITAE_n','cost'});
    best.table = sortrows(best.table, 'cost', 'ascend');
end