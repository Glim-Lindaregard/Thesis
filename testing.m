% AMS_onefile_test_and_export.m
% One-file script that uses your existing:
%   - config()
%   - buildAMS_row(cfg)
%   - visualizeSliderHealth(health, ax)
%
% And exports a PDF with 4 tiles per page (2x2):
%   Row: [Thruster configuration | corresponding AMS]
%
% IMPORTANT FIX (why your "active-only" AMS looked normal):
%   Active (stuck-on) thrusters are NOT a normal bound box; they are FIXED.
%   The attainable wrench set is the healthy-set translated by A_active*u_stuck.
%   Here we implement that translation explicitly while still using buildAMS_row
%   to generate the variable-part facets.
%
% health encoding:
%   1 = healthy (variable)
%   0 = passive failure (off)
%   2 = active failure (stuck on at u_active_stuck)

clear; clc;

%% --- base config ---
cfg_base = config();
A0 = cfg_base.A;
N  = size(A0,2);
assert(N==8,'Expected 8 thrusters.');

u_min_base = cfg_base.u_min(:);
u_max_base = cfg_base.u_max(:);

% buildAMS_row needs this field
cfg_base.N_thrusters = N;

% active stuck level
u_active_stuck = 0.7;

% case limits (same idea as AMS_number_test)
max_failed_total = 4;
max_active       = 3;

% tolerances for optional "zero in hull" classification (copied from AMS_number_test)
tol_inside   = 1e-8;
tol_interior = 1e-6;

%% --- builder ---
if exist('buildAMS_row','file') == 2
    AMS_builder = @buildAMS_row;
elseif exist('buildAMSRow','file') == 2
    AMS_builder = @buildAMSRow;
else
    error('Could not find buildAMS_row or buildAMSRow on MATLAB path.');
end

%% --- enumerate cases (from AMS_number_test, simplified storage) ---
cases = struct('active_idx',{},'passive_idx',{},'health',{},'zero_contains',{},'zero_class',{},'margin',{});

thrusters = 1:N;

count_total = 0; count_zero = 0; count_edge = 0; count_int = 0;

for a = 0:max_active
    active_sets = nchoosek(thrusters, a);
    if a == 0, active_sets = zeros(1,0); end

    for ia = 1:size(active_sets,1)
        active_idx = active_sets(ia,:);
        remaining1 = setdiff(thrusters, active_idx);

        max_passive_here = max_failed_total - a;
        for p = 0:max_passive_here
            passive_sets = nchoosek(remaining1, p);
            if p == 0, passive_sets = zeros(1,0); end

            for ip = 1:size(passive_sets,1)
                passive_idx = passive_sets(ip,:);

                % Build AMS for this case (variable part + translation for active)
                [V, health] = wrenchSet_for_case(cfg_base, AMS_builder, A0, ...
                    u_min_base, u_max_base, u_active_stuck, active_idx, passive_idx);

                % Classify whether zero is in the wrench-space hull (optional but useful)
                [inside, interior, margin] = zeroInsideHull(V, tol_inside, tol_interior);

                if ~inside
                    zclass = "NO";
                elseif interior
                    zclass = "INTERIOR";
                else
                    zclass = "EDGE";
                end

                count_total = count_total + 1;
                if inside
                    count_zero = count_zero + 1;
                    if interior, count_int = count_int + 1; else, count_edge = count_edge + 1; end
                end

                cases(end+1).active_idx    = active_idx; %#ok<SAGROW>
                cases(end).passive_idx     = passive_idx;
                cases(end).health          = health;
                cases(end).zero_contains   = inside;
                cases(end).zero_class      = zclass;
                cases(end).margin          = margin;
            end
        end
    end
end

fprintf('Total cases (<=%d failed, <=%d active): %d', max_failed_total, max_active, count_total);
fprintf('Cases containing zero wrench: %d', count_zero);
fprintf('  - Interior (wrench-hull): %d', count_int);
fprintf('  - Edge (wrench-hull)    : %d', count_edge);

%% --- PDF export: 2 cases per page in 2x2 tiles ---
pdfPath = 'health_ams.pdf';
if exist(pdfPath,'file'), delete(pdfPath); end

% You can filter here if you want (e.g., only interior zero cases):
idx_show = 1:numel(cases);
% idx_show = find(strcmp({cases.zero_class},'INTERIOR')); % uncomment to filter

fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1800 950]);
tlo = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

ci = 1;
while ci <= numel(idx_show)

    for kk = 1:4
        ax = nexttile(tlo, kk);
        cla(ax);
    end

    for r = 1:2
        if ci > numel(idx_show), break; end

        k = idx_show(ci);
        h = cases(k).health;
        act = cases(k).active_idx;
        pas = cases(k).passive_idx;

        % Left: thruster health using your existing function
        axL = nexttile(tlo, 2*r - 1);
        cla(axL);
        visualizeSliderHealth(h, axL);

        title(axL, { ...
            sprintf('case #%d | A=[%s] P=[%s]', k, num2str(act), num2str(pas)), ...
            sprintf('%s | margin=%.3e', cases(k).zero_class, cases(k).margin) ...
            }, 'Interpreter','none');

        % Right: AMS in wrench space (convex hull)
        axR = nexttile(tlo, 2*r);
        cla(axR);

        [V, ~] = wrenchSet_for_case(cfg_base, AMS_builder, A0, ...
            u_min_base, u_max_base, u_active_stuck, act, pas);

        plotHull3_local(axR, V);
        title(axR, 'Attainable Moment Set');
        xlabel(axR,'F_x'); ylabel(axR,'F_y'); zlabel(axR,'T_z');
        grid(axR,'on'); axis(axR,'tight'); view(axR,[55 25]);

        ci = ci + 1;
    end

    drawnow;
    exportgraphics(fig, pdfPath, 'Append', true);
end

close(fig);
fprintf('Done. PDF saved: %s', pdfPath);

%% ===================== helper functions =====================

function [V, health] = wrenchSet_for_case(cfg_base, AMS_builder, A0, umin0, umax0, u_stuck, active_idx, passive_idx)
% Build wrench-space point cloud V (3xM) for a failure case.
% Uses buildAMS_row for the variable part, then applies the correct active translation.

N = size(A0,2);

% health encoding
health = ones(1,N);
health(passive_idx) = 0;
health(active_idx)  = 2;

% --- VARIABLE PART CONFIG ---
% Passive thrusters: removed (A col zero, bounds 0)
% Active thrusters : removed from variable part (A col zero, bounds 0)
% Healthy thrusters: nominal bounds

cfg = cfg_base;
cfg.N_thrusters = N;

lb = umin0; ub = umax0;

if ~isempty(passive_idx)
    lb(passive_idx) = 0; ub(passive_idx) = 0;
end
if ~isempty(active_idx)
    lb(active_idx) = 0; ub(active_idx) = 0;
end

Avar = A0;
Avar(:, passive_idx) = 0;
Avar(:, active_idx)  = 0;

cfg.A     = Avar;
cfg.u_min = lb;
cfg.u_max = ub;

% build AMS facets in thruster space (variable part)
Uraw = AMS_builder(cfg);
Uvar = coerceU(Uraw, N);           % N x M

% --- INSERT FIXED VALUES ---
U = Uvar;
if ~isempty(passive_idx)
    U(passive_idx,:) = 0;
end
if ~isempty(active_idx)
    U(active_idx,:)  = u_stuck;
end

% Map with TRUE A (including active columns)
V = A0 * U;

% de-duplicate to make hull more stable
V = unique(round(V.', 12), 'rows').';
end

function U = coerceU(Uraw, N)
% From AMS_number_test: flatten N x 4 x F into N x (4F)
if isnumeric(Uraw)
    U = Uraw;
    if ndims(U) == 3
        if size(U,1) ~= N
            error('3D AMS points: expected first dim %d, got %d.', N, size(U,1));
        end
        U = reshape(U, N, []);
    end
elseif isstruct(Uraw)
    if isfield(Uraw,'U'), U = Uraw.U;
    elseif isfield(Uraw,'verts'), U = Uraw.verts;
    elseif isfield(Uraw,'points'), U = Uraw.points;
    else
        error('AMS builder returned a struct, but no recognized field (U/verts/points).');
    end
    if isnumeric(U) && ndims(U) == 3
        if size(U,1) ~= N
            error('3D AMS points in struct: expected first dim %d, got %d.', N, size(U,1));
        end
        U = reshape(U, N, []);
    end
else
    error('AMS builder returned unsupported type: %s', class(Uraw));
end

if size(U,1) ~= N && size(U,2) == N
    U = U.';
end
if size(U,1) ~= N
    error('AMS points have wrong dimension. Expected %dxM, got %dx%d.', N, size(U,1), size(U,2));
end

bad = any(~isfinite(U),1);
U(:,bad) = [];
if isempty(U), error('AMS point set is empty after cleaning.'); end
end

function plotHull3_local(ax, V)
cla(ax); hold(ax,'on');
X = V.'; % Mx3

if size(X,1) < 4
    plot3(ax, V(1,:), V(2,:), V(3,:), 'k.');
    hold(ax,'off');
    return;
end

try
    K = convhulln(X);
    trisurf(K, X(:,1), X(:,2), X(:,3), 'Parent', ax, 'FaceAlpha', 0.25, 'EdgeAlpha', 0.35);
catch
    plot3(ax, V(1,:), V(2,:), V(3,:), 'k.');
end

plot3(ax, 0,0,0, 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
hold(ax,'off');
end

function [inside, interior, margin] = zeroInsideHull(V, tol_inside, tol_interior)
% From AMS_number_test
inside = false; interior = false; margin = 0;
M = size(V,2);
if M < 4, return; end
X = V.';
try
    K = convhulln(X);
catch
    return;
end
c = mean(X,1).';
margin = inf;
inside = true;

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
