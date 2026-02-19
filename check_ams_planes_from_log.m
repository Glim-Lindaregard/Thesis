function check_ams_planes_from_log(matFile, ad, adLogCsv)
%CHECK_AMS_PLANES_FROM_LOG  Diagnose AMS plane issues from logged .mat file.
%
% Usage examples:
%   check_ams_planes_from_log("planes_mode0.mat", [1.4 1.4 0.39]);
%   check_ams_planes_from_log("planes_mode0.mat", [], "ad_log.csv");
%
% Expected variables in matFile:
%   - n : (M,3) plane normals
%   - b : (M,1) or (1,M) plane offsets
% Optional (strongly recommended for real "wrong-facing" proof):
%   - V : (Nv,3) vertices of the AMS polytope used to generate planes
%
% Plane convention assumed:
%   n*u <= b  (each row of n is a plane normal)
%
% adLogCsv format (optional):
%   Columns named Fx,Fy,Tau OR three numeric columns (Fx Fy Tau)

    if nargin < 2, ad = []; end
    if nargin < 3, adLogCsv = ""; end

    S = load(matFile);
    assert(isfield(S,"n") && isfield(S,"b"), "MAT file must contain variables 'n' and 'b'.");

    n = double(S.n);
    b = double(S.b(:));
    M = size(n,1);
    assert(size(n,2) == 3, "n must be Mx3.");
    assert(numel(b) == M, "b must have length M.");

    planes = [n, b];

    fprintf("Loaded %d planes from %s\n", M, matFile);

    % ------------------------------
    % 1) Check single a_d (if given)
    % ------------------------------
    if ~isempty(ad)
        ad = double(ad(:)).';
        assert(numel(ad)==3, "ad must be [Fx Fy Tau].");
        vals = n*ad(:) - b;                 % <= 0 is feasible
        [mx, idx] = max(vals);
        mx = max(0, mx);

        fprintf("\nSingle-point check:\n");
        fprintf("  ad = [%.4f %.4f %.4f]\n", ad(1), ad(2), ad(3));
        fprintf("  max_violation = %.6e at plane %d\n", mx, idx);
        fprintf("  worst plane = [%.6f %.6f %.6f  %.6f]\n", planes(idx,1), planes(idx,2), planes(idx,3), planes(idx,4));
    end

    % ---------------------------------------------
    % 2) Check a log of a_d points (if CSV provided)
    % ---------------------------------------------
    if strlength(adLogCsv) > 0
        T = readtable(adLogCsv);

        if all(ismember(["Fx","Fy","Tau"], T.Properties.VariableNames))
            U = [T.Fx, T.Fy, T.Tau];
        else
            U = table2array(T);
            assert(size(U,2) >= 3, "CSV must contain Fx,Fy,Tau columns or at least 3 numeric columns.");
            U = U(:,1:3);
        end

        vals = n*U.' - b;     % (M, N)
        maxPerSample = max(vals, [], 1);
        maxPerSample = max(0, maxPerSample);

        fprintf("\nLog check (%d samples from %s):\n", size(U,1), adLogCsv);
        fprintf("  max violation over all samples: %.6e\n", max(maxPerSample));
        fprintf("  fraction violating: %.2f %%\n", 100*mean(maxPerSample > 1e-9));

        [worstV, j] = max(maxPerSample);
        if worstV > 0
            uWorst = U(j,:);
            valsWorst = n*uWorst(:) - b;
            [mx, idx] = max(valsWorst);
            fprintf("  worst sample index %d: ad=[%.4f %.4f %.4f], max=%.6e at plane %d\n", ...
                j, uWorst(1), uWorst(2), uWorst(3), mx, idx);
            fprintf("  worst plane = [%.6f %.6f %.6f  %.6f]\n", planes(idx,1), planes(idx,2), planes(idx,3), planes(idx,4));
        end
    end

    % ---------------------------------------------------------
    % 3) Heuristic diagnostics on the plane set itself
    % ---------------------------------------------------------
    fprintf("\nPlane diagnostics:\n");

    % 3a) Feasibility of the half-space intersection (polytope)
    % Solve: find any u such that n*u <= b
    % Use linprog with 0 objective.
    f = [0;0;0];
    A = n;
    rhs = b;

    opts = optimoptions("linprog","Display","none");
    [u_feas,~,exitflag] = linprog(f, A, rhs, [], [], [], [], opts);

    if exitflag == 1
        fprintf("  Feasibility: OK (found u=[%.4f %.4f %.4f])\n", u_feas(1), u_feas(2), u_feas(3));
    else
        fprintf("  Feasibility: FAIL (linprog exitflag=%d). Plane set inconsistent!\n", exitflag);
        fprintf("  -> If inconsistent, IPOPT will return infeasible iterates and you'll see violations.\n");
    end

    % 3b) Duplicate / near-duplicate planes (same normal + same b)
    tolN = 1e-6; tolB = 1e-6;
    [dupPairs] = find_near_duplicate_planes(planes, tolN, tolB);
    fprintf("  Near-duplicates (same n, same b): %d pairs\n", size(dupPairs,1));
    if ~isempty(dupPairs)
        disp("  Example duplicate pair (i,j):");
        disp(dupPairs(1,:));
    end

    % 3c) Opposite normals pairs (n_i ≈ -n_j)
    cosTol = 0.9999;  % tighten/loosen as needed
    bTol   = 1e-6;
    oppPairs = find_opposite_normal_pairs(planes, cosTol, bTol);
    fprintf("  Opposite-normal suspicious pairs (n_i≈-n_j and b_i≈b_j): %d pairs\n", size(oppPairs,1));
    if ~isempty(oppPairs)
        i = oppPairs(1,1); j = oppPairs(1,2);
        fprintf("  Example pair i=%d, j=%d\n", i, j);
        fprintf("    plane i = [%.6f %.6f %.6f  %.6f]\n", planes(i,1), planes(i,2), planes(i,3), planes(i,4));
        fprintf("    plane j = [%.6f %.6f %.6f  %.6f]\n", planes(j,1), planes(j,2), planes(j,3), planes(j,4));
        fprintf("  Note: This is a heuristic red flag, not a proof.\n");
    end

    % ---------------------------------------------------------
    % 4) REAL wrong-facing test (requires vertices V)
    % ---------------------------------------------------------
    if isfield(S,"V")
        V = double(S.V);
        assert(size(V,2)==3, "V must be Nv x 3.");
        dots = n*V.';                 % (M, Nv)
        maxDot = max(dots, [], 2);    % (M,1)
        viol = maxDot - b;            % should be <= 0 for ALL planes
        bad = find(viol > 1e-6);

        fprintf("\nVertex-consistency test (using V in mat file):\n");
        fprintf("  bad planes (max(n*v) > b): %d / %d\n", numel(bad), M);

        if ~isempty(bad)
            [worstV, k] = max(viol(bad));
            i = bad(k);
            fprintf("  WORST plane idx=%d, violation=%.6e\n", i, viol(i));
            fprintf("  plane = [%.6f %.6f %.6f  %.6f]\n", planes(i,1), planes(i,2), planes(i,3), planes(i,4));
            fprintf("  max(n*v)=%.6f, b=%.6f\n", maxDot(i), b(i));
            fprintf("  -> This is strong evidence the plane is wrong-facing or built from wrong facet.\n");
        else
            fprintf("  All planes are consistent with the provided vertices.\n");
        end
    else
        fprintf("\nVertex-consistency test skipped: no variable 'V' in the mat file.\n");
        fprintf("  -> If you can log V (vertices) once per mode, this test becomes definitive.\n");
    end
end


function dupPairs = find_near_duplicate_planes(planes, tolN, tolB)
    M = size(planes,1);
    dupPairs = [];
    for i = 1:M
        for j = i+1:M
            if norm(planes(i,1:3) - planes(j,1:3)) < tolN && abs(planes(i,4) - planes(j,4)) < tolB
                dupPairs(end+1,:) = [i,j]; %#ok<AGROW>
            end
        end
    end
end

function oppPairs = find_opposite_normal_pairs(planes, cosTol, bTol)
    n = planes(:,1:3);
    b = planes(:,4);

    nn = n ./ max(vecnorm(n,2,2), 1e-12);

    M = size(planes,1);
    oppPairs = [];
    for i = 1:M
        for j = i+1:M
            c = dot(nn(i,:), nn(j,:));
            if c < -cosTol && abs(b(i) - b(j)) < bTol
                oppPairs(end+1,:) = [i,j]; %#ok<AGROW>
            end
        end
    end
end