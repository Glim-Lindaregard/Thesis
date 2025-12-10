function investigate_active_failures_active()
% INVESTIGATE_ACTIVE_FAILURE_EQUILIBRIA
%   For each active-failure pattern, analyze whether a zero net wrench
%   (Fx=Fy=Tau=0) is achievable with healthy thrusters compensating
%   the stuck-ON ones, and how "hard" the healthy thrusters must work.
%
%   A      : 3 x N thruster effectiveness matrix
%   u_max  : N x 1 vector of max thrust per thruster [N]
%
% State encoding (conceptual):
%   1 = healthy (0 <= u_i <= u_max(i))
%   2 = active-failed (u_i = u_max(i))
%   0 = dead (u_i = 0)  [not used here, but supported in analyze_config]
%
% Usage example:
%   investigate_active_failure_equilibria(A, u_max);


    A = [-1,     0,     1,     0,     1,     0,    -1,     0;
          0,     1,     0,     1,     0,    -1,     0,    -1;
      -0.14,  0.14,  0.14, -0.14, -0.14,  0.14,  0.14, -0.14];


    u_max = 0.7*ones(8,1);

    [m, N] = size(A);
    if m ~= 3
        error('A must be 3xN (Fx,Fy,Tau rows).');
    end
    u_max = u_max(:);
    if numel(u_max) ~= N
        error('u_max must have length N.');
    end

    fprintf('================ ACTIVE FAILURE EQUILIBRIA ==================\n');
    fprintf('N thrusters: %d\n', N);

    % ------ Baseline: no active failures ------
    base_state = ones(N,1);  % all healthy
    fprintf('\nBaseline (no active failures):\n');
    analyze_config(A, u_max, base_state, []);

    % ------ For K = 1..N active failures ------
    for K = 1:N
        fprintf('\n=== %d active failure(s) ===\n', K);
        combos = nchoosek(1:N, K);

        for row = 1:size(combos,1)
            active_failed = combos(row,:);

            state = ones(N,1);      % start as healthy
            state(active_failed) = 2;  % mark these as active-failed

            analyze_config(A, u_max, state, active_failed);
        end
    end
end


function analyze_config(A, u_max, state, active_failed)
% ANALYZE_CONFIG
%   A, u_max as above.
%   state: N x 1, in {0,1,2}
%   active_failed: indices of thrusters that are active-failed (for printing).

    N = numel(state);
    healthy_idx = find(state == 1);
    active_idx  = find(state == 2);
    dead_idx    = find(state == 0); %#ok<NASGU>

    A_H = A(:, healthy_idx);
    A_F = A(:, active_idx);
    uF  = u_max(active_idx);

    % Target for healthy thrusters: A_H * u_H = -A_F * uF
    b = -A_F * uF;

    fprintf('Active-failed thrusters: %s\n', mat2str(active_failed));

    % ----- Case 1: no healthy thrusters -----
    if isempty(healthy_idx)
        residual = norm(A_F * uF, 2);
        if residual < 1e-6
            fprintf('  No healthy thrusters needed, already zero wrench.\n');
        else
            fprintf('  No healthy thrusters. Zero wrench IMPOSSIBLE.\n');
            fprintf('  Residual ||A_F u_F|| = %0.4g\n', residual);
        end
        return;
    end

    % ----- Unconstrained least-squares solution -----
    % Solve A_H u_H ≈ b in least-squares sense
    % (if A_H is full row rank, this is the minimum-norm solution).
    uH_ls = A_H \ b;  %  uses backslash; equivalent to least-squares
    residual_ls = norm(A_H * uH_ls - b, 2);

    % Check bounds for this solution
    uH_min = zeros(size(uH_ls));
    uH_max = u_max(healthy_idx);
    violates_lower = (uH_ls < uH_min - 1e-9);
    violates_upper = (uH_ls > uH_max + 1e-9);

    within_bounds = ~any(violates_lower | violates_upper);

    % Saturation measure: fraction of u_max used by each healthy thruster
    usage = uH_ls ./ uH_max;
    usage(isnan(usage)) = 0; % in case u_max=0 somewhere
    max_usage = max(usage);
    idx_max = healthy_idx(1 + (find(usage == max_usage, 1, 'first')-1));  % map back to global index

    fprintf('  LS residual ||A_H u_H + A_F u_F|| = %0.4g\n', residual_ls);

    if residual_ls > 1e-5
        fprintf('  -> Even unconstrained solution cannot make net wrench zero.\n');
    end

    if within_bounds
        fprintf('  -> Zero wrench is ACHIEVABLE within 0 <= u_H <= u_max.\n');
        fprintf('     Max usage among healthy thrusters: %0.3f (thruster %d)\n', max_usage, idx_max);
        if max_usage > 0.9
            fprintf('     !!! Requires ~full thrust on some healthy thruster -> fragile equilibrium.\n');
        end
    else
        fprintf('  -> LS solution violates bounds:\n');
        if any(violates_lower)
            fprintf('     Below 0 for thrusters: %s\n', mat2str(healthy_idx(violates_lower)));
        end
        if any(violates_upper)
            fprintf('     Above u_max for thrusters: %s\n', mat2str(healthy_idx(violates_upper)));
        end

        % Clip to [0,u_max] and see best achievable residual
        uH_clip = min(max(uH_ls, uH_min), uH_max);
        residual_clip = norm(A_H * uH_clip - b, 2);
        fprintf('     After clipping, residual = %0.4g\n', residual_clip);

        if residual_clip < 1e-3
            fprintf('     -> Approximate zero wrench possible, but at bound(s).\n');
        else
            fprintf('     -> Significant unavoidable wrench -> drift.\n');
        end
    end
end
