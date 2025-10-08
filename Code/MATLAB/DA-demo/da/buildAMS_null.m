function AMS = buildAMS_null(A, umin, umax, tol)
% buildAMS_null - Construct Attainable Moment Set using null space method
%
% Implements the null space-based vertex finding from Tang et al. (2011)
% combined with row-space facet construction (Section 3.2, parts 2 and 4)
%
% Inputs:
%   A     - Control effectiveness matrix (3 x m)
%   umin  - Minimum control limits (m x 1 or 1 x m)
%   umax  - Maximum control limits (m x 1 or 1 x m)
%   tol   - Tolerance for numerical comparisons (optional, default 1e-10)
%
% Output:
%   AMS - Structure containing:
%         .vertices - Cell array of vertex moments (3x1 each)
%         .controls - Cell array of corresponding control vectors
%         .facets   - Structure array with fields:
%                     .Tau    - 4x3 matrix of facet vertices [A;B;C;D]
%                     .normal - 3x1 outward normal vector
%                     .pair   - [i j] indices of varying actuators
%                     .type   - 'max' or 'min'

    if nargin < 4
        tol = 1e-10;
    end
    
    % Ensure row vectors for consistency
    umin = umin(:)';
    umax = umax(:)';
    
    [~, m] = size(A);
    
    fprintf('Building AMS using null-space vertex method...\n');
    
    % Step 1: Find vertices using null space method
    [u_vertices, moment_vertices] = findVerticesNullSpace(A, umin, umax, tol);
    
    n_vert = length(moment_vertices);
    fprintf('  Found %d vertices\n', n_vert);
    
    % Step 2: Construct facets using row-space method (like buildAMS.m)
    % This is the key: use row-space for facets even with null-space vertices
    facets = constructFacetsRowSpace(A, umin, umax, tol);
    
    fprintf('  Built %d facets\n', length(facets));
    
    % Package output
    AMS.vertices = moment_vertices;
    AMS.controls = u_vertices;
    AMS.facets = facets;
    AMS.m = m;
    AMS.num_facets_built = length(facets);
    AMS.num_facets_expected = m * (m - 1);
end

function [u_vertices, moment_vertices] = findVerticesNullSpace(A, umin, umax, tol)
    % Find AMS vertices using null space method (Section 3.2, part 4)
    m = size(A, 2);
    
    % Step 1: Generate all 2^m saturated control vectors
    u_all = generateAllSaturated(umin, umax);
    
    % Step 2: Find null space basis
    N = null(A);  % Null space basis
    
    if isempty(N)
        % No null space - all vertices map uniquely
        u_vertices_mat = u_all;
    else
        % Step 3: Remove vertices that map to interior using Theorem 1
        u_vertices_mat = filterInteriorVertices(A, N, u_all, umin, umax, tol);
    end
    
    % Convert to cell arrays
    n_vert = size(u_vertices_mat, 2);
    u_vertices = cell(n_vert, 1);
    moment_vertices = cell(n_vert, 1);
    
    for i = 1:n_vert
        u_vertices{i} = u_vertices_mat(:, i);
        moment_vertices{i} = A * u_vertices_mat(:, i);
    end
end

function u_all = generateAllSaturated(umin, umax)
    % Generate all 2^m vertices where each control is at min or max
    m = length(umin);
    n_vertices = 2^m;
    u_all = zeros(m, n_vertices);
    
    for i = 0:(n_vertices-1)
        binary = dec2bin(i, m) - '0';  % Convert to binary array
        for j = 1:m
            if binary(j) == 0
                u_all(j, i+1) = umin(j);
            else
                u_all(j, i+1) = umax(j);
            end
        end
    end
end

function u_keep = filterInteriorVertices(A, N, u_all, umin, umax, tol)
    % Remove vertices that map to interior of AMS using null space (Eq. 6, 15, 16)
    m = size(A, 2);
    n_all = size(u_all, 2);
    
    to_remove = false(n_all, 1);
    
    % Method: Use 3-combination approach (more thorough, from paper)
    indices = 1:m;
    combos = nchoosek(indices, 3);
    
    for c = 1:size(combos, 1)
        idx = combos(c, :);
        A_sub = A(:, idx);
        
        % Check if invertible
        if abs(det(A_sub)) < tol
            continue;
        end
        
        % Remaining indices
        idx_rem = setdiff(indices, idx);
        if isempty(idx_rem)
            continue;
        end
        
        A_rem = A(:, idx_rem);
        
        % For each basis direction in remaining dimensions
        for k = 1:length(idx_rem)
            x_k = zeros(length(idx_rem), 1);
            x_k(k) = 1;
            
            % Solve for the 3 selected controls (Eq. 14)
            xi_sub = -A_sub \ (A_rem * x_k);
            
            % Construct full xi vector (null space direction)
            xi = zeros(m, 1);
            xi(idx) = xi_sub;
            xi(idx_rem) = x_k;
            
            % Build vertices that would map to interior (Eq. 15 and 16)
            u_interior_pos = constructNullSpaceVertex(xi, umin, umax, 1);
            u_interior_neg = constructNullSpaceVertex(xi, umin, umax, -1);
            
            % Mark for removal
            for i = 1:n_all
                if norm(u_all(:, i) - u_interior_pos) < tol || ...
                   norm(u_all(:, i) - u_interior_neg) < tol
                    to_remove(i) = true;
                end
            end
        end
    end
    
    % Keep only boundary vertices
    u_keep = u_all(:, ~to_remove);
end

function u_k = constructNullSpaceVertex(xi, umin, umax, sign_choice)
    % Construct vertex according to Eq. (6), (15), or (16)
    % sign_choice: 1 for Eq.(15), -1 for Eq.(16)
    m = length(xi);
    u_k = zeros(m, 1);
    
    for i = 1:m
        if sign_choice * xi(i) > 0
            u_k(i) = umax(i);
        elseif sign_choice * xi(i) < 0
            u_k(i) = umin(i);
        else  % xi(i) == 0
            u_k(i) = umax(i);  % Arbitrary choice when on boundary
        end
    end
end

function facets = constructFacetsRowSpace(A, umin, umax, tol)
    % Construct facets using row-space method (Section 3.2, part 2)
    % This is equivalent to the approach in buildAMS.m
    
    [~, m] = size(A);
    facets = [];
    fcount = 0;
    
    % Iterate all unordered pairs (i,j), i<j
    for i = 1:(m-1)
        ci = A(:, i);
        for j = (i+1):m
            cj = A(:, j);
            
            % Compute facet normal
            n = cross(ci, cj);
            if norm(n) < tol
                continue;  % Skip parallel columns
            end
            
            % Triple-product signs: s_k = sign(n · c_k)
            n_dot_A = n' * A;
            s = sign(n_dot_A);
            s([i j]) = 0;  % Free variables
            
            % Two parallel facets: 'max' (+n) and 'min' (-n)
            for which = 1:2
                if which == 1
                    ftype = 'max';
                    normal = n;
                    choose = +1;
                else
                    ftype = 'min';
                    normal = -n;
                    choose = -1;
                end
                
                % Build template u by locking all k ≠ {i,j}
                u_tmpl = zeros(1, m);
                for k = 1:m
                    if k == i || k == j, continue; end
                    if choose * s(k) > 0
                        u_tmpl(k) = umax(k);
                    else
                        u_tmpl(k) = umin(k);
                    end
                end
                
                % Four vertices: A(min,min), B(max,min), C(max,max), D(min,max)
                Uverts = zeros(4, m);
                Uverts(1, :) = u_tmpl; Uverts(1, i) = umin(i); Uverts(1, j) = umin(j);
                Uverts(2, :) = u_tmpl; Uverts(2, i) = umax(i); Uverts(2, j) = umin(j);
                Uverts(3, :) = u_tmpl; Uverts(3, i) = umax(i); Uverts(3, j) = umax(j);
                Uverts(4, :) = u_tmpl; Uverts(4, i) = umin(i); Uverts(4, j) = umax(j);
                
                % Map to AMS vertices τ = A·u
                Tau = (A * Uverts')';  % 4×3
                
                % Store facet
                fcount = fcount + 1;
                facets(fcount).Tau = Tau;
                facets(fcount).normal = normal(:);
                facets(fcount).pair = [i j];
                facets(fcount).type = ftype;
                facets(fcount).Uverts = Uverts;
            end
        end
    end
end