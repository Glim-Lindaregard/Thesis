function [Uuniq, Uctrl] = d4_failure_classes_and_controllable()
% d4_failure_classes_and_controllable
% -------------------------------------------------------------------------
% Builds D4-unique (rotation + mirroring) 3-state failure classes for 8
% thrusters, then filters the subset whose AMS (via buildAMS_row) is:
%   (i)  inflated (full-dimensional in R^3, i.e., not a sheet)
%   (ii) contains the origin (0,0,0) in wrench space
%
% State encoding per thruster (row vector length 8):
%   0 = passive failed (stuck OFF)    => u_min=0,   u_max=0
%   1 = neutral/healthy              => u_min=0,   u_max=umax
%   2 = active failed (stuck ON)     => u_min=umax,u_max=umax
%
% OUTPUTS
%   Uuniq : Nx8 canonical representatives under D4
%   Uctrl : Mx8 subset of Uuniq that passes the AMS tests
%
% REQUIREMENTS in your path:
%   - config()
%   - buildAMS_row(cfg)   (now patched to be numerically robust)
%
% IMPORTANT:
%   This code assumes config() provides:
%     cfg.pos   (8x2 positions)
%     cfg.beta  (8x1 or 1x8 thruster angles, radians)
%     cfg.A     (3x8 wrench mapping)
%     cfg.u_max (8x1 or 1x8, with nominal umax in entry 1)
% -------------------------------------------------------------------------

    cfg0 = config();

    % Robust D4 permutations from (pos + direction), not pos-only.
    perms = d4_permutations_from_geometry(cfg0.pos, cfg0.beta);

    % Enumerate all 3^8 patterns and keep one canonical representative per D4 class
    Uuniq = d4_unique_3state(perms);

    % Filter to controllable via buildAMS_row(cfg)
    Uctrl = filterControllableClasses_buildAMS(Uuniq);
end

% =====================================================================
% 1) Enumerate + D4-canonicalize
% =====================================================================

function Uuniq = d4_unique_3state(perms)
    Nthr = 8;
    nAll = 3^Nthr;

    seen = containers.Map('KeyType','char','ValueType','logical');
    Uuniq = zeros(nAll, Nthr);
    nU = 0;

    for idx = 0:(nAll-1)
        s = base3_vec(idx, Nthr);               % 0..2
        sCanon = canonical_under_d4(s, perms);  % 0..2

        key = vec_key(sCanon);
        if ~isKey(seen, key)
            seen(key) = true;
            nU = nU + 1;
            Uuniq(nU,:) = sCanon;
        end
    end

    Uuniq = Uuniq(1:nU,:);
end

function sCanon = canonical_under_d4(s, perms)
% Choose canonical representative as minimal base-3 integer among all D4 images
    best = s;
    bestVal = base3_int(s);

    for g = 1:size(perms,1)
        sg = apply_perm(s, perms(g,:));
        v = base3_int(sg);
        if v < bestVal
            bestVal = v;
            best = sg;
        end
    end
    sCanon = best;
end

function sg = apply_perm(s, perm)
% perm(i)=j means element i moves to j (destination index)
    sg = zeros(size(s));
    for i = 1:numel(s)
        sg(perm(i)) = s(i);
    end
end

function v = base3_int(s)
% s is 0..2, interpret as base-3 number (s(1) least significant)
    v = 0; p = 1;
    for i = 1:numel(s)
        v = v + s(i) * p;
        p = p * 3;
    end
end

function s = base3_vec(idx, n)
% idx -> length-n base-3 digits in {0,1,2}, least significant first
    s = zeros(1,n);
    for i = 1:n
        s(i) = mod(idx, 3);
        idx = floor(idx/3);
    end
end

function key = vec_key(s)
% stable string key for containers.Map
    key = char(s + '0'); % '0','1','2'
end

% =====================================================================
% 2) D4 permutations from geometry (pos + direction)
% =====================================================================

function perms = d4_permutations_from_geometry(pos, beta)
% Build 8 D4 permutations by transforming both:
%   - position pos(i,:) in body frame
%   - force direction f(i,:) = [cos(beta) sin(beta)]
%
% This avoids ambiguity if two thrusters share the same corner position.

    pos = pos(:,1:2);
    beta = beta(:);
    f = [cos(beta) sin(beta)];
    f = f ./ vecnorm(f,2,2);

    T = d4_transforms(); % 1x8 cell of 2x2
    tol = 1e-10;

    N = size(pos,1);
    perms = zeros(8, N);
    row = 0;

    for k = 1:numel(T)
        Tk = T{k};
        pos2 = (Tk * pos')';
        f2   = (Tk * f')';

        perm = zeros(1,N);
        used = false(1,N);

        for i = 1:N
            % Match transformed thruster i to some original j
            posErr = vecnorm(pos - pos2(i,:), 2, 2);
            dirErr = vecnorm(f   - f2(i,:),   2, 2);
            score  = posErr + dirErr;

            % pick best unused match
            [~, order] = sort(score, 'ascend');
            j = order(1);
            if used(j) && numel(order) > 1
                j = order(find(~used(order), 1, 'first'));
            end

            if isempty(j) || score(j) > tol
                error('D4 mapping failed: no geometry match within tol. Check cfg.pos/cfg.beta ordering.');
            end

            perm(i) = j;
            used(j) = true;
        end

        if numel(unique(perm)) ~= N
            error('D4 mapping not bijective. Geometry matching ambiguous; check cfg.pos/cfg.beta.');
        end

        row = row + 1;
        perms(row,:) = perm;
    end
end

function T = d4_transforms()
% 8 D4 transforms as 2x2 matrices acting on [x;y]
    R0   = [ 1  0;  0  1];
    R90  = [ 0 -1;  1  0];
    R180 = [-1  0;  0 -1];
    R270 = [ 0  1; -1  0];

    % reflection about x-axis, then rotate
    F0   = [ 1  0;  0 -1];
    F90  = R90  * F0;
    F180 = R180 * F0;
    F270 = R270 * F0;

    T = {R0, R90, R180, R270, F0, F90, F180, F270};
end

% =====================================================================
% 3) Controllability filter using buildAMS_row(cfg)
% =====================================================================

function Uctrl = filterControllableClasses_buildAMS(U)
% For each class s in U:
%   - set cfg.u_min/u_max per your 0/1/2 definition
%   - Facets = buildAMS_row(cfg)
%   - Map facet vertices to wrench vertices (A * Ufacet)
%   - Check inflated (3D volume) and contains origin

    cfg0 = config();
    umax_nom = cfg0.u_max(1);

    % tolerances (scale-aware)
    tol_round = 1e-9;
    tol_rank  = 1e-12;
    tol_vol   = 1e-12;   % will be scaled per case
    tol_in    = 1e-9;    % will be scaled per case

    Uctrl = zeros(size(U));
    m = 0;

    for kcase = 1:size(U,1)
        s = U(kcase,:);

        cfg = cfg0;
        umin = zeros(8,1);
        umax = zeros(8,1);

        for i = 1:8
            if s(i) == 0
                umin(i) = 0.0;
                umax(i) = 0.0;
            elseif s(i) == 1
                umin(i) = 0.0;
                umax(i) = umax_nom;
            elseif s(i) == 2
                umin(i) = umax_nom;
                umax(i) = umax_nom;
            else
                error('Invalid state (must be 0/1/2).');
            end
        end

        cfg.u_min = umin;
        cfg.u_max = umax;

        % --- Build facets in thruster space
        Facets = buildAMS_row(cfg); % 8 x 4 x nf

        % --- Convert to wrench-space vertices
        W = facets_to_wrench_vertices(cfg.A, Facets); % M x 3

        if isempty(W) || size(W,1) < 4 || any(~isfinite(W(:)))
            continue
        end

        % de-dup / stabilize
        W = unique(round(W / tol_round) * tol_round, 'rows');

        if size(W,1) < 4
            continue
        end

        % full-dimensional check
        if rank(W - mean(W,1), tol_rank) < 3
            continue
        end

        % scale-aware tolerances
        scale = max(1, max(abs(W), [], 'all'));
        tolV  = max(tol_vol, 1e-12 * scale^3);
        tolH  = max(tol_in,  1e-9  * scale);

        % convex hull volume
        try
            [K, vol] = convhulln(W);
        catch
            continue
        end
        r = origin_inradius(W, K);          % margin (meters*newtons-ish)
        if r <= 1e-6 * max(1, max(abs(W),[],'all'))   % tune
            continue
        end


        % Require origin to be strictly inside (not on boundary)
        if ~origin_in_hull_interior(W, K)
            continue
        end


        m = m + 1;
        Uctrl(m,:) = s;
    end

    Uctrl = Uctrl(1:m,:);
end

function W = facets_to_wrench_vertices(A, Facets)
% Facets: 8 x 4 x nf  (each facet gives 4 vertices in thruster space)
% Returns W: (4*nf) x 3 wrench vertices
    nf = size(Facets, 3);
    if nf == 0
        W = [];
        return
    end

    W = zeros(4*nf, 3);
    t = 0;

    for k = 1:nf
        Uk = Facets(:,:,k);     % 8x4 vertices (thruster space)
        Vk = (A * Uk).';        % 4x3 vertices (wrench space)
        W(t+1:t+4,:) = Vk;
        t = t + 4;
    end
end

function inside = origin_in_hull_halfspace(W, K, tol)
% Half-space test for origin using hull triangles K.
% Uses c = mean(W) as guaranteed interior point (convex combination).
    p = [0 0 0];
    c = mean(W,1);

    for f = 1:size(K,1)
        v0 = W(K(f,1),:);
        v1 = W(K(f,2),:);
        v2 = W(K(f,3),:);

        n = cross(v1 - v0, v2 - v0);
        nn = norm(n);
        if nn < tol
            continue
        end

        % orient outward: outward normal points away from interior point c
        if dot(n, c - v0) > 0
            n = -n;
        end

        % inside if dot(n, p - v0) <= tol for all outward facets
        if dot(n, p - v0) > tol
            inside = false;
            return
        end
    end

    inside = true;
end
function inside = origin_in_hull_interior(W, K)
% Strict interior test for p = [0 0 0] in conv(W)
% Returns true only if origin is at least eps away from every hull face.

    p = [0 0 0];
    c = mean(W,1);                % guaranteed inside conv(W)
    scale = max(1, max(abs(W), [], 'all'));
    eps_in = 1e-6 * scale;        % tune if you want stricter/looser

    for f = 1:size(K,1)
        v0 = W(K(f,1),:);
        v1 = W(K(f,2),:);
        v2 = W(K(f,3),:);

        n = cross(v1 - v0, v2 - v0);
        nn = norm(n);
        if nn < 1e-12 * scale
            continue
        end

        % orient outward using interior point c
        if dot(n, c - v0) > 0
            n = -n;
        end

        % Strict inside: origin must be strictly on the inside side
        % i.e. dot(n, p - v0) <= -eps
        if dot(n, p - v0) > -eps_in
            inside = false;
            return
        end
    end

    inside = true;
end
function r = origin_inradius(W, K)
% Returns the inradius at the origin for the convex hull conv(W).
% r > 0  => origin strictly inside
% r = 0  => origin on boundary (or outside due to numeric issues)
%
% Uses hull triangles K and orients each face outward using an interior point.

    p = [0 0 0];
    c = mean(W,1);  % convex combo of vertices => inside conv(W)
    scale = max(1, max(abs(W), [], 'all'));
    epsn = 1e-14 * scale;

    r = inf;

    for f = 1:size(K,1)
        v0 = W(K(f,1),:);
        v1 = W(K(f,2),:);
        v2 = W(K(f,3),:);

        n = cross(v1 - v0, v2 - v0);
        nn = norm(n);
        if nn < epsn, continue, end

        % orient outward
        if dot(n, c - v0) > 0
            n = -n;
        end

        % Face inequality: n·x <= n·v0  (for outward n)
        b = dot(n, v0);

        % distance from origin to plane along outward normal is b/||n||
        d = b / nn;

        r = min(r, d);
        if r <= 0
            r = 0;
            return
        end
    end

    if ~isfinite(r), r = 0; end
end
