function exportZeroClassHealthPDF(cases, uniqKeys, g, zeroIdx, pdfPath, ...
    cfg_base, AMS_builder, u_min_base, u_max_base, u_active_stuck, N)
% Export 2 cases per page in a 2x2 grid:
%   Row 1: [Health | AMS]
%   Row 2: [Health | AMS]
%
% Filters to INTERIOR only (recommended for clean PDF).

if nargin < 5 || isempty(pdfPath)
    pdfPath = 'zero_classes.pdf';
end
if exist(pdfPath,'file')
    delete(pdfPath);
end

% Representative for each symmetry class
numClassesAll = numel(uniqKeys);
repsAll = zeros(numClassesAll,1);
for ci = 1:numClassesAll
    members = zeroIdx(g == ci);
    repsAll(ci) = members(1);
end

% ---- Filter to INTERIOR only ----
isInterior = ([cases(repsAll).zero_class] == "INTERIOR");
reps = repsAll(isInterior);

fprintf('Exporting %d / %d symmetry-class reps (INTERIOR only) to %s\n', ...
    numel(reps), numel(repsAll), pdfPath);

% quick preview toggle
exportOnlyFirstPage = false;  % set true for 1 page (2 cases)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exportOnlyFirstPage
    reps = reps(1:min(2,numel(reps)));
end

if isempty(reps)
    warning('No INTERIOR representatives to export.');
    return;
end

% Stable figure size => stable PDF pages
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1800 950]);
tlo = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

ci = 1;
numCases = numel(reps);
pageNr = 1;
while ci <= numCases

    % clear page tiles
    for kk = 1:4
        ax = nexttile(tlo, kk);
        cla(ax);
    end

    % 2 rows => 2 cases/page
    for r = 1:2
        if ci > numCases, break; end

        idx = reps(ci);
        h   = cases(idx).health;

        % Left tile: health plot
        axL = nexttile(tlo, 2*r - 1);
        cla(axL);
        visualizeSliderHealth(h, axL);

        act = cases(idx).active_idx;
        pas = cases(idx).passive_idx;

        title(axL, { ...
            sprintf('rep case #%d | A=[%s] P=[%s]', idx, num2str(act), num2str(pas)), ...
            sprintf('INTERIOR | margin=%.3e', cases(idx).margin) ...
            }, 'Interpreter','none', 'Color',[0 0.6 0]);

        % Right tile: AMS plot
        axR = nexttile(tlo, 2*r);
        cla(axR);

        % --- NEW "buildAMS" for our purposes ---
        [Ucurrent, cfg_case] = buildAMS_case(cfg_base, AMS_builder, ...
            u_min_base, u_max_base, u_active_stuck, act, pas, N);

        % --- NEW visualizeAMS that draws into axR (no new figure) ---
        visualizeAMS_ax(Ucurrent, cfg_case.Acurrent, "normal", 0.7*[1 1 1]', [], [], axR);

        ci = ci + 1;
        print("Page number %d done", pageNr);
        pageNr = pageNr +1;
    end

    drawnow;
    exportgraphics(fig, pdfPath, 'Append', true);
end

close(fig);
fprintf('Done. PDF saved: %s\n', pdfPath);
end


function [Ucurrent, cfg_case] = buildAMS_case(cfg_base, AMS_builder, ...
    u_min_base, u_max_base, u_active_stuck, active_idx, passive_idx, N)


% Default: nominal bounds (healthy thrusters behave normally)
lb = u_min_base(:);
ub = u_max_base(:);

% Passive failures: clamp OFF
if ~isempty(passive_idx)
    lb(passive_idx) = 0.0;
    ub(passive_idx) = 0.0;
end

% Active failures: clamp ON (stuck-on)
if ~isempty(active_idx)
    lb(active_idx) = u_active_stuck;
    ub(active_idx) = u_active_stuck;
end

% END FIX ----------------------------------------------------


% Passive failures are already OFF by default

% END FIX ----------------------------------------------------

cfg_case = cfg_base;
cfg_case.u_min = lb;
cfg_case.u_max = ub;

if ~isfield(cfg_case,'Acurrent') || isempty(cfg_case.Acurrent)
    cfg_case.Acurrent = cfg_case.A;
end

Ucurrent = AMS_builder(cfg_case);

% --- FORCE 8x4xF format for visualizeAMS_ax ---
if isnumeric(Ucurrent) && ndims(Ucurrent) == 2
    % if builder returned flattened 8x(4F), reshape to 8x4xF
    if mod(size(Ucurrent,2), 4) ~= 0
        error('AMS_builder returned %dx%d, not divisible by 4 columns -> cannot interpret facets.', ...
              size(Ucurrent,1), size(Ucurrent,2));
    end
    Ucurrent = reshape(Ucurrent, size(Ucurrent,1), 4, []);
end

end



function visualizeAMS_ax(U, Asys, mode, aProduced, ad, abc, ax)
% visualizeAMS_ax
% Copy of your visualizeAMS, modified to draw into an axes (tile) instead of creating figures.

if nargin < 4 || isempty(aProduced), aProduced = [0;0;0]; end
if nargin < 5 || isempty(ad),        ad        = [0;0;0]; end
if nargin < 6 || isempty(abc),       abc       = [0;0;0]; end
if nargin < 7 || isempty(ax) || ~isgraphics(ax,'axes')
    error('visualizeAMS_ax requires a valid axes handle as 7th argument.');
end

mode = lower(string(mode));

% --- preset styles (same as your function) ---
switch mode
    case "export"
        FaceColor       = [0.8 0.95 0.95];
        EdgeColor       = [0.3 0.3 0.3];
        FaceAlpha       = 0.95;
        EdgeAlpha       = 1;
        LineWidth       = 0.8;
        BackgroundColor = [1 1 1];
        GridColor       = 0.3*[1 1 1];
        GridAlpha       = 1;
        GridStyle       = '--';
        FontSize        = 18;
        ViewAngles      = [55 30];
        ShowProduced    = true;
        ShowDesired     = true;
        ShowBasis       = true;
    case "normal"
        FaceColor       = 0.75*[0.8 0.95 0.95];
        EdgeColor       = [0 0 0];
        FaceAlpha       = 0.95;
        EdgeAlpha       = 1;
        LineWidth       = 0.8;
        BackgroundColor = [1 1 1];          % <- IMPORTANT: white for PDF tiles
        GridColor       = 0.2*[1 1 1];
        GridAlpha       = 1;
        GridStyle       = '--';
        FontSize        = 12;
        ViewAngles      = [55 25];
        ShowProduced    = true;
        ShowDesired     = true;
        ShowBasis       = false;
    otherwise
        error('visualizeAMS_ax:UnknownMode','Mode must be "normal" or "export".');
end

facetIndex = 1;

% --- axes setup (NO FIGURE CREATION) ---
cla(ax);
hold(ax,'on');
set(ancestor(ax,'figure'),'Renderer','painters','Color',BackgroundColor);
set(ax,'Color',BackgroundColor);

count  = size(U,3);
scale  = 2;
center = [0 0 0];

axis(ax,'manual'); axis(ax,'equal'); axis(ax,'vis3d');
xlim(ax,[-scale,scale]);
ylim(ax,[-scale,scale]);
zlim(ax,[-scale/4,scale/4]);
view(ax, ViewAngles(1), ViewAngles(2));

lgHandles = [];
lgLabels  = {};

% --- draw facets ---
for k = 1:count
    Uk    = U(:,:,k);
    Vk    = Asys*Uk;
    verts = Vk;

    A = verts(:,1)'; B = verts(:,2)'; C = verts(:,3)'; D = verts(:,4)';
    tris = [A; B; C; D];
    f    = [1 2 3 4];

    isHighlighted = (k == facetIndex) && (ShowBasis || ShowProduced);

    if isHighlighted
        patch(ax,'Faces',f,'Vertices',tris, ...
          'FaceColor','r','FaceAlpha',0.5, ...
          'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
          'EdgeLighting','none','FaceLighting','flat','HandleVisibility','off');
    else
        patch(ax,'Faces',f,'Vertices',tris, ...
          'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
          'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
          'EdgeLighting','none','FaceLighting','flat','HandleVisibility','off');
    end
end

% --- produced / desired wrenches ---
if ShowProduced
    hProd = quiver3(ax, center(1),center(2),center(3), ...
                    aProduced(1), aProduced(2), aProduced(3), ...
                    'AutoScale','off','Color','b', ...
                    'LineWidth',1.5,'MaxHeadSize',0.8);
    lgHandles(end+1) = hProd;
    lgLabels{end+1}  = 'Moment produced';
end

if ShowDesired
    hDes = quiver3(ax, center(1),center(2),center(3), ...
                   ad(1), ad(2), ad(3), ...
                   'AutoScale','off','Color','r', ...
                   'LineWidth',2,'MaxHeadSize',0.8);
    lgHandles(end+1) = hDes;
    lgLabels{end+1}  = 'Desired moment';
end

% --- axes / labels / legend ---
box(ax,'on'); grid(ax,'on');
ax.GridColor     = GridColor;
ax.GridAlpha     = GridAlpha;
ax.GridLineStyle = GridStyle;
ax.XMinorGrid    = 'off';
ax.YMinorGrid    = 'off';
ax.ZMinorGrid    = 'off';
ax.FontSize      = FontSize;
ax.LineWidth     = LineWidth;
ax.Projection    = 'perspective';
ax.XColor        = 'k';
ax.YColor        = 'k';
ax.ZColor        = 'k';

xlabel(ax,'F_x (N)','Interpreter','tex');
ylabel(ax,'F_y (N)','Interpreter','tex');
zlabel(ax,'T_z (N·m)','Interpreter','tex');
title(ax,'Attainable Moment Set','Color','k');

if ~isempty(lgHandles)
    legend(ax, lgHandles, lgLabels{:}, 'TextColor','k', 'Location','best');
end

hold(ax,'off');
end
