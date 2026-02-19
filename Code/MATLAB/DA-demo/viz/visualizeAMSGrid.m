function visualizeAMSGrid(uCache, failureTime, opts, aProduced, ad)
% VISUALIZEAMSGRID  Grid of AMS plots over failure sequence defined by failureTime.
% Inputs
%   uCache       : 1x256 struct, fields uCache(k).U (m×4×F), uCache(k).A (3×m)
%                  Indexing uses alive-bit mask: LSB = thruster 1. idx = 1+bin2dec(alivebits).
%   failureTime  : 1×8 double. Time [s] each thruster fails. Inf => never; t<0 => failed at start.
%   opts         : (optional) struct with fields similar to your visualizeAMS (colors, etc.)
%   aProduced    : (optional) 3×1 vector, drawn in every panel if provided
%   ad           : (optional) 3×1 vector, drawn in every panel if provided
%
% Notes
%   - Builds scenarios at t = [-Inf, unique(finite failure times in ascending order)].
%   - In scenario j, thruster i is alive iff t < failureTime(i).
%   - All subplots share same axis limits and view.
%
% Glim: no camera spin, no index highlighting, uses patch.

    if nargin < 3 || isempty(opts), opts = struct(); end

    FaceColor       = getfield_with_default(opts,'FaceColor',[0.78 0.92 0.92]);
    EdgeColor       = getfield_with_default(opts,'EdgeColor','k');
    FaceAlpha       = getfield_with_default(opts,'FaceAlpha',0.95);
    EdgeAlpha       = getfield_with_default(opts,'EdgeAlpha',0.85);
    LineWidth       = getfield_with_default(opts,'LineWidth',0.75);
    BackgroundColor = getfield_with_default(opts,'BackgroundColor','w');
    GridColor       = getfield_with_default(opts,'GridColor',[0.55 0.55 0.55]);
    GridAlpha       = getfield_with_default(opts,'GridAlpha',0.35);
    GridStyle       = getfield_with_default(opts,'GridLineStyle','--');
    FontSize        = getfield_with_default(opts,'FontSize',11);
    ViewAzEl        = getfield_with_default(opts,'View',[45 30]);

    % --- Build scenario thresholds: start + each finite failure time ---
    ft = double(failureTime(:)).';
    ft(~isfinite(ft)) = inf;    % NaN -> never
    thresholds = [-inf, unique(ft(isfinite(ft)))];   % ascending

    % --- For each scenario, compute alive bits -> uCache index ---
    S = numel(thresholds);                  % <= 9
    idx = zeros(1,S);
    aliveBits = false(S,8);                 % for subplot titles/debug

    for s = 1:S
        tcut = thresholds(s);
        alive = (tcut < ft);                % alive if it hasn't failed yet
        aliveBits(s,:) = alive;
        idx(s) = 1 + sum(uint16(alive).*uint16(2.^(0:7)));  % LSB=thruster1
    end

    % --- Precompute global axes from all scenarios for consistency ---
    allXYZ = [];   % collect all vertices across panels
    panels = cell(1,S);  % hold panel data (U,A) to avoid double indexing later
    Tcache4 = ...
    [
    1 0 0 0 0 1 1 1;
    0 1 0 0 1 1 0 1;
    1 0 1 1 0 1 0 0;
    0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0;
    1 0 0 1 0 1 0 0;
    0 0 1 0 0 1 1 1;
    1 1 0 0 1 1 0 0;
    0 1 0 1 0 1 0 1;
    0 0 0 1 0 1 1 1;
    1 0 0 1 0 1 1 0;
    0 0 1 1 0 1 0 1;
    0 0 0 0 1 1 1 1;
    0 1 1 0 0 1 1 0;
    ];

    for s = 1:S
        k = Tcache4(s,:);
        entry = uCache(k);
        panels{s}.U = entry.U;   % m×4×F
        panels{s}.A = entry.A;   % 3×m

        U = panels{s}.U; A = panels{s}.A;
        F = size(U,3);
        for f = 1:F
            Vk = A * U(:,:,f);                   % 3×4
            allXYZ = [allXYZ; Vk.'];             %#ok<AGROW>  % 4×3
        end
    end

    if isempty(allXYZ), allXYZ = [0 0 0; 1 0 0; 0 1 0; 0 0 1]; end
    [xmin,xmax] = bounds(allXYZ(:,1));
    [ymin,ymax] = bounds(allXYZ(:,2));
    [zmin,zmax] = bounds(allXYZ(:,3));

    % Pad a bit for arrows
    pad = 0.01 * max([abs(xmax-xmin), abs(ymax-ymin), abs(zmax-zmin), 1]);
    xl = [xmin-pad, xmax+pad];
    yl = [ymin-pad, ymax+pad];
    zl = [zmin-pad, zmax+pad];

    % --- Make a near-square grid ---
    [nRows, nCols] = bestTile(S);

    % --- Figure ---
    fig = figure('Color',BackgroundColor, 'Name','AMS Grid'); %#ok<NASGU>
    set(gcf, 'Position', [80 80 1920 1080]);  % make the window bigger


    for s = 1:S
        subplot(nRows, nCols, s); hold on; box on; grid on;
        ax = gca;  %#ok<NASGU>
        ax.GridColor = GridColor; ax.GridAlpha = GridAlpha; ax.GridLineStyle = GridStyle;
        ax.MinorGridLineStyle = GridStyle; ax.XMinorGrid='on'; ax.YMinorGrid='on'; ax.ZMinorGrid='on';
        ax.FontSize = FontSize; axis equal; axis vis3d;
        xlim(xl); ylim(yl); zlim(zl);
        view(ViewAzEl(1), ViewAzEl(2));

        U = panels{s}.U; A = panels{s}.A; F = size(U,3);

        % Facets
        for f = 1:F
            Vk = A * U(:,:,f);   % 3×4
            verts = Vk.';        % 4×3 (A,B,C,D rows)
            patch('Faces',[1 2 3 4], 'Vertices', verts, ...
                  'FaceColor',FaceColor, 'FaceAlpha',FaceAlpha, ...
                  'EdgeColor',EdgeColor, 'EdgeAlpha',EdgeAlpha, ...
                  'LineWidth',LineWidth, 'EdgeLighting','none', ...
                  'FaceLighting','flat', 'HandleVisibility','off');
        end

        % Optional arrows (draw last so they're on top)
        ctr = [0 0 0];
        if nargin >= 4 && ~isempty(aProduced)
            quiver3(ctr(1),ctr(2),ctr(3), aProduced(1),aProduced(2),aProduced(3), ...
                'AutoScale','off','Color',[0.2 0.35 0.9],'LineWidth',1.4,'MaxHeadSize',0.8);
        end
        if nargin >= 5 && ~isempty(ad)
            quiver3(ctr(1),ctr(2),ctr(3), ad(1),ad(2),ad(3), ...
                'AutoScale','off','Color',[0.15 0.65 0.15],'LineWidth',1.2,'MaxHeadSize',0.8);
        end

        xlabel('Fx'); ylabel('Fy'); zlabel('Tau');

        % Title indicates which thrusters are alive (1) / failed (0)
        bits = aliveBits(s,:);
        title(sprintf('Scenario %d  |  alive bits = [%s]', ...
              s, sprintf('%d', bits(end:-1:1))), 'Interpreter','none'); 
        % Note: shown MSB..LSB for readability; logic uses LSB=thruster1
    end
end

% ---------- helpers ----------
function [r,c] = bestTile(n)
% near-square tiling up to 9 panels (works for any n)
    c = ceil(sqrt(n));
    r = ceil(n / c);
end

function v = getfield_with_default(s,name,def)
    if isfield(s,name), v = s.(name); else, v = def; end
end
