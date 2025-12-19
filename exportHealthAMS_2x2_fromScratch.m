% AMS_onefile_test_and_export.m
% One-file script that:
%   1) Defines slider geometry (A-matrix) and thruster layout
%   2) Generates failure cases (passive + active)
%   3) Builds and plots AMS in wrench space with CORRECT active-failure modeling
%      (active stuck-on thrusters = translation of the variable set)
%   4) Exports a PDF with 4 figures per page (2x2):
%        Row i: [Thruster configuration | Corresponding AMS]
%
% This is rebuilt from scratch and does NOT rely on buildAMSRow/buildAMS_row.
% It will work as intended even if your old builder ignored u_min.
%
% HOW TO USE
%   - Put this file on your MATLAB path.
%   - Run it. It writes health_ams.pdf in the current folder.
%
% Notes:
%   health encoding:
%     1 = healthy (variable, 0..u_max)
%     0 = passive failure (off, fixed 0)
%     2 = active failure (stuck on, fixed u_stuck)

clear; clc;

%% ---------------------------- CONFIG ----------------------------
% Platform half-size (for plotting only)
a_plot = 0.22;                 % meters-ish scale (only affects left plot)

% Thruster force magnitude bounds
u_min_base = zeros(8,1);
u_max_base = 0.7*ones(8,1);

% Active (stuck-on) thrust level
u_stuck = 0.7;

% Sampling density per thruster for AMS point-cloud approximation
% (2^8 = 256 points if levels=[0,1] only; more levels => more points)
levels_per_thruster = 2;       % set 2 for fast and robust; 3+ for smoother hull

% Export
pdfPath = 'health_ams.pdf';
if exist(pdfPath,'file'), delete(pdfPath); end

%% ----------------------- THRUSTER GEOMETRY ----------------------
% Define 8 thrusters, two per corner of a square, in body frame.
% Each thruster produces a force along a fixed direction, and contributes
% torque tau = r x F (scalar in 2D).

% Corner positions (square)
r = 0.195;                      % distance to corner-ish (m)
% Use a slightly off-corner placement to resemble your layout
pos = [
     r, -0.140;  % T1
   0.140, -r;    % T2
    -r, -0.140;  % T3
  -0.140, -r;    % T4
    -r,  0.140;  % T5
  -0.140,  r;    % T6
     r,  0.140;  % T7
   0.140,  r;    % T8
];

% Force directions (BODY). These must match your A-matrix convention.
% Pick a common "two per corner" pattern. You can swap/rotate if your real
% A differs; the logic still works.
% Here: thrusters roughly push outward from the body.
dir = [
  1,  0;   % T1 pushes +x
  0, -1;   % T2 pushes -y
 -1,  0;   % T3 pushes -x
  0, -1;   % T4 pushes -y
 -1,  0;   % T5 pushes -x
  0,  1;   % T6 pushes +y
  1,  0;   % T7 pushes +x
  0,  1;   % T8 pushes +y
];

% Build A = [Fx; Fy; tau] mapping, size 3x8
A0 = zeros(3,8);
for i = 1:8
    fx = dir(i,1);
    fy = dir(i,2);
    tau = pos(i,1)*fy - pos(i,2)*fx;   % r x F (z-component)
    A0(:,i) = [fx; fy; tau];
end

%% ---------------------------- CASES -----------------------------
% Rebuild "AMS_number_test" style enumeration: group by number of failures
% and show unique-ish cases. For simplicity, we generate all combinations
% up to max_failed, with optional active/passive modes.

max_failed = 3;                % adjust to taste
include_active = true;
include_passive = true;

cases = build_failure_cases(8, max_failed, include_passive, include_active);

fprintf('Generated %d cases (max_failed=%d, active=%d, passive=%d)\n', ...
    numel(cases), max_failed, include_active, include_passive);

%% ----------------------- PDF 2x2 EXPORT -------------------------
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1800 950]);
tlo = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

ci = 1;
while ci <= numel(cases)

    % clear tiles
    for k = 1:4
        ax = nexttile(tlo,k);
        cla(ax);
    end

    for row = 1:2
        if ci > numel(cases), break; end

        c = cases(ci);
        h = c.health;
        passive_idx = find(h==0);
        active_idx  = find(h==2);

        % Left tile: thruster health
        axL = nexttile(tlo, 2*row-1);
        visualizeSliderHealth_local(axL, a_plot, pos, dir, h);
        title(axL, sprintf('case %d | active=[%s] passive=[%s]', ...
            ci, num2str(active_idx), num2str(passive_idx)), 'Interpreter','none');

        % Right tile: AMS in wrench space (correct active modeling)
        axR = nexttile(tlo, 2*row);
        V = compute_AMS_wrench(A0, u_min_base, u_max_base, h, u_stuck, levels_per_thruster);
        plotHull3_local(axR, V);
        title(axR, 'AMS (wrench space)');
        xlabel(axR,'F_x'); ylabel(axR,'F_y'); zlabel(axR,'\tau_z');
        grid(axR,'on'); axis(axR,'tight'); view(axR,[55 25]);

        ci = ci + 1;
    end

    drawnow;
    exportgraphics(fig, pdfPath, 'Append', true);
end

close(fig);
fprintf('Done. PDF saved: %s\n', pdfPath);

%% ============================ FUNCTIONS ============================

function cases = build_failure_cases(N, max_failed, include_passive, include_active)
% Returns struct array with field .health (1xN)
% health values: 1 healthy, 0 passive, 2 active

cases = struct('health', {});

% Always include baseline
cases(end+1).health = ones(1,N);

for k = 1:max_failed
    combs = nchoosek(1:N, k);
    for ci = 1:size(combs,1)
        idx = combs(ci,:);

        % Passive-only case
        if include_passive
            h = ones(1,N);
            h(idx) = 0;
            cases(end+1).health = h;
        end

        % Active-only case
        if include_active
            h = ones(1,N);
            h(idx) = 2;
            cases(end+1).health = h;
        end

        % Mixed (some active, some passive) for this failure set
        if include_active && include_passive && k>=2
            % simple mixed split: first half active, rest passive
            split = floor(k/2);
            h = ones(1,N);
            h(idx(1:split)) = 2;
            h(idx(split+1:end)) = 0;
            cases(end+1).health = h;
        end
    end
end
end

function V = compute_AMS_wrench(A0, u_min_base, u_max_base, health, u_stuck, levels_per_thruster)
% Compute achievable wrench set V = A0*u, using correct modeling:
%  - healthy thrusters: u in [u_min, u_max] (typically [0, umax])
%  - passive failed:    u fixed 0
%  - active failed:     u fixed u_stuck
%
% The set with fixed components is a translation of the variable set.
% We implement it directly by sampling the healthy dimensions and inserting
% fixed values for failed thrusters.

N = size(A0,2);
health = health(:).';

healthy_idx = find(health==1);
passive_idx = find(health==0);
active_idx  = find(health==2);

% Build per-thruster value grids
vals = cell(1,N);
for i = 1:N
    if any(i==passive_idx)
        vals{i} = 0.0;
    elseif any(i==active_idx)
        vals{i} = u_stuck;
    else
        % healthy
        umin = u_min_base(i);
        umax = u_max_base(i);
        if levels_per_thruster <= 2
            vals{i} = [umin, umax];
        else
            vals{i} = linspace(umin, umax, levels_per_thruster);
        end
    end
end

% Cartesian product of vals{1}..vals{N}
U = cartprod_numeric(vals);  % returns N x M

% Map to wrench space
V = A0 * U;                  % 3 x M

% Remove duplicates (optional but makes hull more stable for small grids)
V = unique(round(V.', 12), 'rows').';
end

function U = cartprod_numeric(vals)
% vals: 1xN cell, each contains 1xMi numeric vector
% output U: N x prod(Mi)

N = numel(vals);
M = 1;
for i=1:N
    M = M * numel(vals{i});
end
U = zeros(N, M);

rep = 1;
block = M;
for i = 1:N
    vi = vals{i}(:).';
    mi = numel(vi);
    block = block / mi;
    pattern = repelem(vi, block);
    fullvec = repmat(pattern, 1, rep);
    U(i,:) = fullvec;
    rep = rep * mi;
end
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
    trisurf(K, X(:,1), X(:,2), X(:,3), ...
        'Parent', ax, ...
        'FaceAlpha', 0.25, ...
        'EdgeAlpha', 0.35);
catch
    plot3(ax, V(1,:), V(2,:), V(3,:), 'k.');
end

plot3(ax, 0,0,0, 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);

hold(ax,'off');
end

function visualizeSliderHealth_local(ax, a, pos, dir, health)
% Draw square body and thrusters, color by health.
% Arrows point INTO the thruster location (direction of generated force).

health = health(:).';
if numel(health) ~= 8, error('health must have 8 elements'); end

cla(ax); hold(ax,'on'); axis(ax,'equal');
set(ax,'Color','w','XColor','k','YColor','k');

xlim(ax, [-1.2*a, 1.2*a]);
ylim(ax, [-1.2*a, 1.2*a]);
rectangle(ax,'Position',[-a -a 2*a 2*a],'LineWidth',2);

colHealthy = [0.35 0.35 0.35];
colPassive = [1.00 0.20 0.20];
colActive  = [0.20 0.45 1.00];

L = 0.25*a;
for i = 1:8
    if health(i) == 1
        col = colHealthy;
    elseif health(i) == 0
        col = colPassive;
    else
        col = colActive;
    end

    plot(ax, pos(i,1), pos(i,2), 'ko', 'MarkerFaceColor', col, 'MarkerSize',10);

    d = dir(i,:);
    d = d / max(norm(d), eps);

    % Arrow should point INTO thruster location along generated force direction
    p_end   = pos(i,:);
    p_start = p_end - L*d;

    quiver(ax, p_start(1), p_start(2), (p_end(1)-p_start(1)), (p_end(2)-p_start(2)), ...
        0, 'Color', col, 'LineWidth',2, 'MaxHeadSize',0.6);
end

axis(ax,'off');
hold(ax,'off');
end
