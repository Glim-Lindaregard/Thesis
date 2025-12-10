% Assume A is already in the workspace, size 3 x m

A = [-1,     0,     1,     0,     1,     0,    -1,     0;
      0,     1,     0,     1,     0,    -1,     0,    -1;
  -0.14,  0.14,  0.14, -0.14, -0.14,  0.14,  0.14, -0.14];


u_max = 0.7;

% === INPUTS ==========================================================
% Assume A (3 x m) and u_max (either scalar or 1 x m / m x 1) exist.
% We assume u_min = 0 for all thrusters (unidirectional thrusters).
% =====================================================================

[rows, m] = size(A);
if rows ~= 3
    error('This script assumes A is 3 x m (Fx; Fy; Tau rows).');
end

% Build u_max vector
if isscalar(u_max)
    u_max_vec = repmat(u_max, 1, m);
else
    u_max_vec = u_max(:).';  % force row vector
    if numel(u_max_vec) ~= m
        error('u_max must be scalar or have length = number of columns of A.');
    end
end

thrusters = 1:m;

% ---------- helper to compute bounds for given A_sub and u_max_sub -----
compute_bounds = @(A_sub, umax_sub) ...
    local_compute_bounds(A_sub, umax_sub);

% ---------- baseline (no failures) ------------------------------------
[Fx0_min, Fx0_max, Fy0_min, Fy0_max, Tau0_min, Tau0_max] = ...
    compute_bounds(A, u_max_vec);

fprintf('Baseline (no failures):\n');
fprintf('  Fx in [%.3f, %.3f]\n', Fx0_min, Fx0_max);
fprintf('  Fy in [%.3f, %.3f]\n', Fy0_min, Fy0_max);
fprintf('  Tau in [%.3f, %.3f]\n\n', Tau0_min, Tau0_max);

% Flags for which directions were available originally
orig_plusFx  = (Fx0_max  > 0);
orig_minusFx = (Fx0_min  < 0);
orig_plusFy  = (Fy0_max  > 0);
orig_minusFy = (Fy0_min  < 0);
orig_plusTau = (Tau0_max > 0);
orig_minusTau= (Tau0_min < 0);

% ---------- sweep over number of failed thrusters ---------------------
for r = 1:m
    combFail = nchoosek(thrusters, r);
    fprintf('=== Failing %d thruster(s) ===\n', r);
    
    for i = 1:size(combFail, 1)
        failed = combFail(i, :);
        
        % Alive set
        aliveMask = true(1, m);
        aliveMask(failed) = false;
        alive = find(aliveMask);
        
        A_alive    = A(:, alive);
        umax_alive = u_max_vec(alive);
        
        % Rank of alive effectiveness matrix
        rA = rank(A_alive);
        
        % Bounds with this failure pattern
        [Fx_min, Fx_max, Fy_min, Fy_max, Tau_min, Tau_max] = ...
            compute_bounds(A_alive, umax_alive);
        
        % Which directions are still possible?
        now_plusFx   = (Fx_max  > 0);
        now_minusFx  = (Fx_min  < 0);
        now_plusFy   = (Fy_max  > 0);
        now_minusFy  = (Fy_min  < 0);
        now_plusTau  = (Tau_max > 0);
        now_minusTau = (Tau_min < 0);
        
        % Lost directions compared to baseline
        lost_dirs = {};
        if orig_plusFx && ~now_plusFx
            lost_dirs{end+1} = '+Fx'; %#ok<*AGROW>
        end
        if orig_minusFx && ~now_minusFx
            lost_dirs{end+1} = '-Fx';
        end
        if orig_plusFy && ~now_plusFy
            lost_dirs{end+1} = '+Fy';
        end
        if orig_minusFy && ~now_minusFy
            lost_dirs{end+1} = '-Fy';
        end
        if orig_plusTau && ~now_plusTau
            lost_dirs{end+1} = '+Tau';
        end
        if orig_minusTau && ~now_minusTau
            lost_dirs{end+1} = '-Tau';
        end
        
        if isempty(lost_dirs)
            lost_str = 'none';
        else
            lost_str = strjoin(lost_dirs, ' ');
        end
        
        fprintf('  Failed: %s | Alive: %s | rank(A_alive) = %d | lost: %s\n', ...
            mat2str(failed), mat2str(alive), rA, lost_str);
    end
    
    fprintf('\n');
end

% ================= local function =====================================
function [Fx_min, Fx_max, Fy_min, Fy_max, Tau_min, Tau_max] = ...
    local_compute_bounds(A_sub, umax_sub)
% A_sub : 3 x k
% umax_sub : 1 x k, u in [0, umax_sub]

    if isempty(A_sub)
        Fx_min = 0; Fx_max = 0;
        Fy_min = 0; Fy_max = 0;
        Tau_min = 0; Tau_max = 0;
        return;
    end

    % ensure row
    umax = umax_sub(:).';
    
    % Row for Fx, Fy, Tau
    rowFx  = A_sub(1, :);
    rowFy  = A_sub(2, :);
    rowTau = A_sub(3, :);
    
    [Fx_min, Fx_max]   = one_row_bounds(rowFx,  umax);
    [Fy_min, Fy_max]   = one_row_bounds(rowFy,  umax);
    [Tau_min, Tau_max] = one_row_bounds(rowTau, umax);
end

function [amin, amax] = one_row_bounds(a_row, umax)
% For u in [0, umax], compute min/max of sum_j a_j * u_j.
% u_min = 0 assumed.

    a = a_row(:).';
    amin = 0;
    amax = 0;
    for j = 1:numel(a)
        aj = a(j);
        ujmax = umax(j);
        if aj > 0
            % contribution in [0, aj*ujmax]
            amax = amax + aj * ujmax;
        elseif aj < 0
            % contribution in [aj*ujmax, 0]
            amin = amin + aj * ujmax;
        end
    end
end
