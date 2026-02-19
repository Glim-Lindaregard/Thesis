%% healthy_planes_from_8x4x56_facets.m
% MATLAB clone of Python simplified AMS cache builder (healthy only)
% facets: 8x4x56 (Uk per facet), A: 3x8

clear; clc;

% ---- REQUIREMENTS in workspace (or load them) ----
% A      : 3x8
% facets : 8x4x56
%
% Example:
% load('healthy_ams.mat','A','facets');

c = config();

A = c.A;

facets = buildAMS_row(c);

assert(exist('A','var')==1, 'A not found in workspace.');
assert(exist('facets','var')==1, 'facets not found in workspace.');
assert(isequal(size(A), [3 8]), 'A must be 3x8.');
assert(ndims(facets)==3 && size(facets,1)==8 && size(facets,2)==4, 'facets must be 8x4xK.');

K = size(facets,3);

% ---- Tolerances (match Python defaults) ----
tol_n = 1e-6;     % direction tolerance: (1 - dot) <= tol_n
tol_b = 1e-6;     % b update tolerance
eps_n = 1e-12;    % degenerate normal threshold

%% --------- Build vertex cloud + per-facet normals ----------
V_all = zeros(4*K, 3);    % each facet contributes 4 vertices in wrench space
normals = zeros(K, 3);
valid = false(K,1);

row = 1;
for k = 1:K
    Uk = facets(:,:,k);      % 8x4
    Vk = A * Uk;             % 3x4  (wrench-space quad)

    % store vertices (4x3) into V_all
    V_all(row:row+3, :) = Vk.';   % (4x3)
    row = row + 4;

    % robust normal from quad
    try
        n = robust_normal_from_quad(Vk, eps_n);  % 1x3 unit
        normals(k,:) = n;
        valid(k) = true;
    catch
        % degenerate facet -> ignore its normal, but keep vertices in V_all
        valid(k) = false;
    end
end

% Keep all vertices; keep only valid normals
V_all = V_all(1:row-1, :);          % (Nv x 3)
normals = normals(valid, :);        % (K_valid x 3)

if isempty(normals)
    error('No valid facet normals computed.');
end

% interior-ish point (convex combo of vertices)
c = mean(V_all, 1);                 % 1x3

%% --------- Cluster normals by signed direction + compute supporting b ----------
cluster_n = zeros(0,3);
cluster_b = zeros(0,1);

for i = 1:size(normals,1)
    n = normals(i,:);               % 1x3

    % Supporting b for this direction (use all vertices)
    vals = V_all * n.';             % (Nv x 1)
    b = max(vals);                  % scalar

    % Orient outward: ensure interior point satisfies n*c' <= b
    [n, b] = orient_outward(n, b, c);

    % Assign to an existing cluster if same signed direction
    assigned = false;
    for j = 1:size(cluster_n,1)
        if same_direction_signed(cluster_n(j,:), n, tol_n)
            % same signed direction -> keep tightest support (max b)
            if b > cluster_b(j) + tol_b
                cluster_b(j) = b;
                cluster_n(j,:) = n;
            end
            assigned = true;
            break;
        end
    end

    % New cluster
    if ~assigned
        cluster_n(end+1,:) = n; %#ok<AGROW>
        cluster_b(end+1,1) = b; %#ok<AGROW>
    end
end

% Outputs (match your Python naming)
Nmat = cluster_n;   % (Np x 3)
bvec = cluster_b;   % (Np x 1)

fprintf('Healthy planes found: %d\n', size(Nmat,1));

%% --------- Definitive sanity checks (containment + outwardness) ----------
epsContain = 1e-7;
margin = 1e-9;

viol = 0;
for k = 1:size(Nmat,1)
    n = Nmat(k,:);
    b = bvec(k);

    worst = max(V_all * n.' - b);     % should be <= 0
    inside = dot(n, c) - b;           % should be < 0

    if worst > epsContain
        fprintf('Containment FAIL plane %d: max(n^T v - b)=%.3e\n', k, worst);
        viol = viol + 1;
    end
    if inside > -margin
        fprintf('Outward FAIL plane %d: n^T c - b=%.3e\n', k, inside);
        viol = viol + 1;
    end
end

if viol == 0
    fprintf('All planes passed containment + outward checks.\n');
else
    fprintf('Total violations: %d\n', viol);
end


%% --------- Plot facets + attached normals ---------

figure; hold on; grid on; axis equal;
xlabel('Fx'); ylabel('Fy'); zlabel('\tau');
title('Healthy AMS Facets with Attached Normals');
view(35,20);

% Plot full vertex cloud (light gray)
scatter3(V_all(:,1), V_all(:,2), V_all(:,3), 8, [0.7 0.7 0.7], 'filled');

scale = 0.15 * max(vecnorm(V_all - mean(V_all),2,2));  % arrow length scale

for k = 1:K
    Uk = facets(:,:,k);        % 8x4
    Vk = A * Uk;               % 3x4
    Vfacet = Vk.';             % 4x3

    % Compute facet normal
    try
        n = robust_normal_from_quad(Vk, eps_n);
    catch
        continue
    end

    % Supporting b and orientation (same logic as before)
    vals = V_all * n.';
    [b, jmax] = max(vals);
    x0 = V_all(jmax,:);   % support point on plane

    [n, b] = orient_outward(n, b, mean(V_all,1));

    % --- Plot facet ---
    patch('Vertices', Vfacet, ...
          'Faces', [1 2 3 4], ...
          'FaceAlpha', 0.15, ...
          'EdgeColor', 'k');

    % --- Plot attached normal ---
    quiver3(x0(1), x0(2), x0(3), ...
            scale*n(1), scale*n(2), scale*n(3), ...
            0, 'LineWidth', 1.5, 'Color', 'r');
end
%% --------- OPTIONAL: quick plot of vertices ----------
figure; hold on; grid on; axis equal;
scatter3(V_all(:,1), V_all(:,2), V_all(:,3), 10, 'filled');
plot3(c(1), c(2), c(3), 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
xlabel('Fx'); ylabel('Fy'); zlabel('\tau');
title('Healthy AMS vertex cloud and interior point');
view(35,20);

%% ================= LOCAL FUNCTIONS =================

function n = unit_vec(v, eps)
    nv = norm(v);
    if nv < eps
        error('Degenerate vector (near zero norm).');
    end
    n = v / nv;
end

function n = robust_normal_from_quad(Vk, eps)
    % Vk: 3x4 quad vertices; returns 1x3 unit normal (sign not canonicalized)
    if ~isequal(size(Vk), [3 4])
        error('Expected Vk shape 3x4, got %dx%d', size(Vk,1), size(Vk,2));
    end

    v = Vk.'; % 4x3
    tri = [1 2 3;
           1 3 4;
           1 2 4;
           2 3 4];

    for t = 1:4
        i = tri(t,1); j = tri(t,2); k = tri(t,3);
        n_raw = cross(v(j,:) - v(i,:), v(k,:) - v(i,:));
        if norm(n_raw) > eps
            n = unit_vec(n_raw, eps);
            return;
        end
    end
    error('Degenerate facet: cannot compute stable normal.');
end

function tf = same_direction_signed(n1, n2, tol_dir)
    % Signed direction match: dot ~ +1 (does NOT treat n and -n as same)
    d = dot(n1, n2);
    tf = (1.0 - d) <= tol_dir;
end

function [n_out, b_out] = orient_outward(n, b, c)
    % Ensure interior point satisfies n^T c <= b; otherwise flip both
    if dot(n, c) > b
        n_out = -n;
        b_out = -b;
    else
        n_out = n;
        b_out = b;
    end
end