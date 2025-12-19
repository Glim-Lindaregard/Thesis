%% exportInteriorZeroAMS_symmetry_pdf.m
% Exports AMS (wrench-space hull) to a multipage PDF for INTERIOR-zero cases,
% grouped by D4 symmetry. One representative per symmetry class by default.

clear; clc;

%% --- base config + AMS builder ---
cfg_base = config();
A = cfg_base.A;
N = size(A,2);
assert(N == 8, 'Assumes 8 thrusters.');



u_min_base = cfg_base.u_min(:);
u_max_base = cfg_base.u_max(:);

if exist('buildAMS_row','file') == 2
    AMS_builder = @buildAMS_row;
elseif exist('buildAMSRow','file') == 2
    AMS_builder = @buildAMSRow;
else
    error('Need buildAMS_row or buildAMSRow on path.');
end

% --- nominal stats (no failures) ---
cfg0 = cfg_base;
U0 = AMS_builder(cfg0);
U0 = coerceU(U0, N);
V0 = A * U0;              % 3xM

nom = amsStats(V0);       % struct with vol, Fx, Fy, Tz


%% --- constraints / knobs ---
max_failed_total = 4;
max_active       = 3;
u_active_stuck   = 0.7;

tol_inside   = 1e-10;
tol_interior = 1e-6;

% Output PDF
out_pdf = fullfile(pwd, 'AMS_interior_zero_symmetry.pdf');
if exist(out_pdf,'file'), delete(out_pdf); end

% Export selection:
EXPORT_ONE_REP_PER_CLASS = true;   % false -> export every interior case found

%% --- thruster geometry for symmetry (D4) ---
thr_pos = [ ...
     195, -140;  % 1
     140, -195;  % 2
    -195, -140;  % 3
    -140, -195;  % 4
    -195,  140;  % 5
    -140,  195;  % 6
     195,  140;  % 7
     140,  195]; % 8

thr_beta_deg = [0;270;180;270;180;90;0;90];
perms = buildD4Perms(thr_pos, thr_beta_deg);

%% --- enumerate and collect INTERIOR-zero cases ---
thrusters = 1:N;

cases = struct('active_idx',{},'passive_idx',{},'health',{},'key',{},'margin',{});

for a = 0:max_active
    if a == 0
        active_sets = zeros(1,0);
    else
        active_sets = nchoosek(thrusters, a);
    end

    for ia = 1:size(active_sets,1)
        active_idx = active_sets(ia,:);
        remaining  = setdiff(thrusters, active_idx);

        max_passive_here = max_failed_total - a;

        for p = 0:max_passive_here
            if p == 0
                passive_sets = zeros(1,0);
            else
                passive_sets = nchoosek(remaining, p);
            end

            for ip = 1:size(passive_sets,1)
                passive_idx = passive_sets(ip,:);

                % bounds for this case
                lb = u_min_base;
                ub = u_max_base;

                if ~isempty(passive_idx)
                    lb(passive_idx) = 0; ub(passive_idx) = 0;
                end
                if ~isempty(active_idx)
                    lb(active_idx) = u_active_stuck;
                    ub(active_idx) = u_active_stuck;
                end

                % build AMS points (thruster space)
                cfg = cfg_base;
                cfg.u_min = lb;
                cfg.u_max = ub;

                Uraw = AMS_builder(cfg);
                U = coerceU(Uraw, N);
                V = A * U;


                [isInterior, margin] = zeroInteriorHull(V, tol_inside, tol_interior);
                if ~isInterior
                    continue;
                end

                h = ones(1,N);
                h(passive_idx) = 0;
                h(active_idx)  = 2;

                cases(end+1).active_idx = active_idx; %#ok<SAGROW>
                cases(end).passive_idx  = passive_idx;
                cases(end).health       = h;
                cases(end).key          = canonicalHealthKey(h, perms);
                cases(end).margin       = margin;
            end
        end
    end
end

if isempty(cases)
    fprintf('No INTERIOR-zero cases found.\n');
    return;
end

%% --- choose which cases to export ---
keys = string({cases.key});
[uniqKeys, ~, grp] = unique(keys);

export_idx = [];
if EXPORT_ONE_REP_PER_CLASS
    for g = 1:numel(uniqKeys)
        members = find(grp == g);
        export_idx(end+1) = members(1); %#ok<SAGROW>
    end
else
    export_idx = 1:numel(cases);
end

fprintf('Found %d interior cases in %d symmetry class(es).\n', numel(cases), numel(uniqKeys));
fprintf('Exporting %d AMS plot(s) to:\n  %s\n', numel(export_idx), out_pdf);

%% PATCH: headless PDF export (no visible figures)
% Drop-in replacement for the export loop + plotting.
% This keeps figures off-screen and never pops up windows.

% --- ensure no figure windows are shown ---
oldVis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

% If you want to be extra safe in newer MATLABs:
% set(0,'DefaultFigureWindowStyle','normal');  % avoid docked popups

%% --- export 4 AMS per page (headless) ---
if exist(out_pdf,'file'), delete(out_pdf); end

pageFig = [];
ax = gobjects(4,1);

for k = 1:2%numel(export_idx)

    % start a new page every 2 cases (=> 4 tiles: AMS+health, AMS+health)
    if mod(k-1,2) == 0
        % flush previous page
        if ~isempty(pageFig)
            exportgraphics(pageFig, out_pdf, 'Append', true, ...
                'ContentType','vector', 'Padding', 5);
            close(pageFig);
        end
    
        % new hidden page figure
        pageFig = figure('Visible','off','Color','w','Renderer','painters');
    
        % 2 cases per page -> 2 rows, each row: [AMS | Health]
        tl = tiledlayout(pageFig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
    end


    % pick case
    c = cases(export_idx(k));

    % rebuild bounds
    lb = u_min_base; ub = u_max_base;
    if ~isempty(c.passive_idx)
        lb(c.passive_idx) = 0; ub(c.passive_idx) = 0;
    end
    if ~isempty(c.active_idx)
        lb(c.active_idx) = u_active_stuck; ub(c.active_idx) = u_active_stuck;
    end

    cfg = cfg_base;
    cfg.u_min = lb;
    cfg.u_max = ub;

    % build AMS (IMPORTANT: keep facet cube N x 4 x F)
    Uraw = AMS_builder(cfg);
    Ufacets = Uraw;
    if isstruct(Uraw)
        if isfield(Uraw,'U'),        Ufacets = Uraw.U;
        elseif isfield(Uraw,'verts'),Ufacets = Uraw.verts;
        elseif isfield(Uraw,'points'),Ufacets = Uraw.points;
        end
    end

    pair = mod(k-1,2);          % 0 or 1 within the page
    healthTile = pair*2 + 1;    % 1 or 3  (LEFT)
    amsTile    = healthTile + 1; % 2 or 4 (RIGHT)

    
    axAMS    = nexttile(tl, amsTile);
    axHealth = nexttile(tl, healthTile);

        % --- fixed AMS axis scaling ---
    Fx_lim  = 1.6;   % pick numbers that cover your data nicely
    Fy_lim  = 1.6;
    Tz_lim  = 0.8;

   
    visualizeSliderHealth(c.health, axHealth);
    title(axHealth,'');
    
    visualizeAMS_ax(axAMS, Ufacets, A, "export");

    % stats
    Upts = coerceU(Ufacets, N);
    Vpts = A * Upts;
    s = amsStats(Vpts);
    
    flags = {};
    if s.Fx  < 0.70*nom.Fx,  flags{end+1} = 'Fx weak'; end
    if s.Fy  < 0.70*nom.Fy,  flags{end+1} = 'Fy weak'; end
    if s.Tz  < 0.70*nom.Tz,  flags{end+1} = '\tau weak'; end   % <-- TeX-safe
    if s.vol < 0.70*nom.vol, flags{end+1} = 'Vol low'; end
    if isempty(flags), flags = {'OK'}; end
    
    title(axAMS, sprintf('AMS volume = %.1f%% of nominal', 100*s.vol/nom.vol), ...
        'FontWeight','bold','Color','k');
    
    % finalize axis geometry FIRST
    xlim(axAMS, [-Fx_lim Fx_lim]);
    ylim(axAMS, [-Fy_lim Fy_lim]);
    zlim(axAMS, [-Tz_lim Tz_lim]);
    axis(axAMS,'vis3d');
    daspect(axAMS,[1 1 1]);
    view(axAMS,35,20);
    camva(axAMS, 9);
    
    % draw red axis-cross
    xl = xlim(axAMS); yl = ylim(axAMS); zl = zlim(axAMS);
    hold(axAMS,'on');
    plot3(axAMS, [xl(1) xl(2)], [0 0], [0 0], 'r-', 'LineWidth', 1.5);
    plot3(axAMS, [0 0], [yl(1) yl(2)], [0 0], 'r-', 'LineWidth', 1.5);
    plot3(axAMS, [0 0], [0 0], [zl(1) zl(2)], 'r-', 'LineWidth', 1.5);
    hold(axAMS,'off');
    
    % NOW place the metrics box (uses final OuterPosition)
    flags_txt = strjoin(flags, newline);
    txt = sprintf(['\\bfReach\\rm\n' ...
                   '|Fx|_{max}: %.3g (%.0f%%)\n' ...
                   '|Fy|_{max}: %.3g (%.0f%%)\n' ...
                   '|\\tau|_{max}: %.3g (%.0f%%)\n\n' ...
                   '\\bfFlags\\rm\n%s'], ...
                   s.Fx, 100*s.Fx/nom.Fx, ...
                   s.Fy, 100*s.Fy/nom.Fy, ...
                   s.Tz, 100*s.Tz/nom.Tz, ...
                   flags_txt);
    
    fig = ancestor(axAMS,'figure');
    p = axAMS.Position;   % [x y w h]

    gap  = 0.01;
    wbox = 0.16;
    hbox = 0.55*p(4);
    
    % try put box to the RIGHT of the axes
    x_try = p(1) + p(3) + gap;
    x_max = 0.99 - wbox;
    
    if x_try <= x_max
        x = x_try;                 % fits -> put outside
    else
        % doesn't fit -> place INSIDE top-right corner with a small inset
        inset = 0.01;
        x = p(1) + p(3) - wbox - inset;
    end
    
    % vertical position (top-aligned), then nudge down
    y = p(2) + p(4) - hbox - 0.03;

    tag = sprintf('AMS_METRICS_BOX_%d', pair);
    delete(findall(fig,'Type','textboxshape','Tag',tag));
    
    annotation(fig,'textbox',[x y wbox hbox], ...
        'String',txt,'Units','normalized','Interpreter','tex', ...
        'FontSize',9,'BackgroundColor','w', ...
        'EdgeColor',0.8*[0.2 0.2 0.2], ...
        'Margin',4,'FitBoxToText','on', ...
        'Tag',tag,'Color',0.1*[1 1 1]);





    
    % title(axAMS, sprintf('key=%s | act=[%s] pas=[%s]', ...
    %     c.key, num2str(c.active_idx), num2str(c.passive_idx)), ...
    %     'Interpreter','none','Color',[0 0 0]);
        
    
    fprintf("case %d done\n",k)

end

% flush last page
if ~isempty(pageFig)
    exportgraphics(pageFig, out_pdf, 'Append', true, ...
        'ContentType','vector', 'Padding', 5);
    close(pageFig);
end

fprintf('Done. Wrote: %s\n', out_pdf);


% restore default visibility
set(0,'DefaultFigureVisible', oldVis);

fprintf('Headless AMS PDF export finished: %s\n', out_pdf);

%% ===================== helpers =====================

function U = coerceU(Uraw, N)
% Return N x M numeric AMS points, flattening N x 4 x F to N x (4F).
    if isnumeric(Uraw)
        U = Uraw;
        if ndims(U) == 3
            if size(U,1) ~= N, error('AMS: expected %d x 4 x F.', N); end
            U = reshape(U, N, []);
        end
    elseif isstruct(Uraw)
        if isfield(Uraw,'U'), U = Uraw.U;
        elseif isfield(Uraw,'verts'), U = Uraw.verts;
        elseif isfield(Uraw,'points'), U = Uraw.points;
        else, error('AMS struct: no U/verts/points field.'); end
        if isnumeric(U) && ndims(U) == 3
            if size(U,1) ~= N, error('AMS: expected %d x 4 x F.', N); end
            U = reshape(U, N, []);
        end
    else
        error('AMS builder returned unsupported type: %s', class(Uraw));
    end

    if size(U,1) ~= N && size(U,2) == N, U = U.'; end
    if size(U,1) ~= N, error('AMS points must be %dxM.', N); end

    bad = any(~isfinite(U),1);
    U(:,bad) = [];
    if isempty(U), error('AMS points empty after cleaning.'); end
end

function [interior, margin] = zeroInteriorHull(V, tol_inside, tol_interior)
% Strict interior test of origin in convex hull of wrench points V (3xM).
    interior = false;
    margin = 0;

    if size(V,2) < 4, return; end
    X = V.'; % Mx3

    try
        K = convhulln(X);
    catch
        return;
    end

    c = mean(X,1).';
    inside = true;
    margin = inf;

    for i = 1:size(K,1)
        p1 = X(K(i,1),:).';
        p2 = X(K(i,2),:).';
        p3 = X(K(i,3),:).';

        n = cross(p2-p1, p3-p1);
        nn = norm(n);
        if nn < 1e-12, continue; end
        n = n/nn;

        if dot(n, c - p1) > 0
            n = -n;
        end

        d = -dot(n, p1);
        if d > tol_inside
            inside = false;
            margin = 0;
            break;
        end
        margin = min(margin, -d);
    end

    interior = inside && (margin > tol_interior);
end

function plotWrenchHull(V)
% Minimal wrench-hull plot (no dependencies on your visualizeAMS).
    X = V.'; % Mx3
    K = convhulln(X);

    trisurf(K, X(:,1), X(:,2), X(:,3), ...
        'FaceAlpha', 0.15, 'EdgeAlpha', 0.25);
    hold on; grid on; axis equal;
    xlabel('F_x'); ylabel('F_y'); zlabel('\tau');

    plot3(0,0,0,'k.','MarkerSize',20);

    % reasonable view
    view(35, 20);
end

function perms = buildD4Perms(pos_xy, beta_deg)
% 8 D4 symmetries: R0,R90,R180,R270 and mirrorX∘Rk
    pos = pos_xy;
    dir = [cosd(beta_deg), sind(beta_deg)];

    rots = cat(3, ...
        [ 1  0; 0  1], ...
        [ 0 -1; 1  0], ...
        [-1  0; 0 -1], ...
        [ 0  1;-1  0]);

    mirrorX = [1 0; 0 -1];

    perms = zeros(8,8);
    g = 0;
    for k = 1:4
        R = rots(:,:,k);
        g = g+1; perms(g,:) = matchPerm(pos, dir, R, pos, dir);
        g = g+1; perms(g,:) = matchPerm(pos, dir, mirrorX*R, pos, dir);
    end
end

function p = matchPerm(pos_src, dir_src, T, pos_tgt, dir_tgt)
    N = size(pos_src,1);
    p = zeros(1,N);
    used = false(1,N);

    tol_pos = 1e-6;
    tol_dir = 1e-6;

    for i = 1:N
        pos_i = (T * pos_src(i,:).').';
        dir_i = (T * dir_src(i,:).').';
        if norm(dir_i) > 0, dir_i = dir_i / norm(dir_i); end

        found = 0;
        for j = 1:N
            if used(j), continue; end
            if norm(pos_tgt(j,:) - pos_i) <= tol_pos
                dj = dir_tgt(j,:);
                if norm(dj) > 0, dj = dj / norm(dj); end
                if norm(dj - dir_i) <= tol_dir || norm(dj + dir_i) <= tol_dir
                    found = j; break;
                end
            end
        end

        if found == 0
            error('D4 match failed at thruster %d. Check geometry/numbering.', i);
        end

        p(i) = found;
        used(found) = true;
    end
end

function key = canonicalHealthKey(h, perms)
% Lexicographically minimal string under D4 permutations.
    best = "";
    for g = 1:size(perms,1)
        ht = applyPerm(h, perms(g,:));
        s = string(sprintf('%d', ht));
        if best == "" || s < best
            best = s;
        end
    end
    key = best;
end

function hnew = applyPerm(h, p)
% p(i)=j means thruster i maps to j
    N = numel(h);
    hnew = zeros(1,N);
    for i = 1:N
        hnew(p(i)) = h(i);
    end
end
function s = amsStats(V)
% V: 3×M wrench points
    X = V.';  % M×3

    % axis reaches
    s.Fx = max(abs(V(1,:)));
    s.Fy = max(abs(V(2,:)));
    s.Tz = max(abs(V(3,:)));

    % volume (convex hull)
    s.vol = 0;
    if size(X,1) >= 4
        try
            [~, s.vol] = convhulln(X);  % returns volume in 3D
        catch
            s.vol = 0;
        end
    end
end
