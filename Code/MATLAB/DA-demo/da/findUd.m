function [ud, found, index, weights] = findUd(AMS, ad, tol)
% findUd  Tang-style facet interpolation (barycentric) to get u from desired wrench ad.
%
% Inputs
%   AMS   : struct with AMS.facets(k).Tau (3x3; rows are facet vertex wrenches)
%                          and AMS.facets(k).Uverts (3xN; rows are vertex control vectors)
%   ad    : 3x1 desired wrench [Fx; Fy; Tau]
%   tol   : (optional) tolerance for inside-triangle check (default 1e-9)
%
% Outputs
%   ud       : Nx1 continuous control (duty fractions/throttles in [0, umax])
%   found    : 1 if a containing facet was found, else 0
%   index    : facet index if found, else 0
%   weights  : [a;b;c] barycentric weights if found, else []

    if nargin < 3 || isempty(tol), tol = 1e-9; end

    found = 0;
    index = 0;
    weights = [];

    % Infer number of thrusters N from first facet
    if isempty(AMS.facets)
        error('AMS.facets is empty.');
    end
    Usample = AMS.facets(1).Uverts;  % expected 3 x N
    N = size(Usample, 2);
    ud = zeros(N,1);

    for k = 1:numel(AMS.facets)
        Tau_k = AMS.facets(k).Tau;     % 3x3, each row is a vertex wrench (1x3)
        U_k   = AMS.facets(k).Uverts;  % 3xN, each row is a vertex control (1xN)

        % Basic shape checks
        if size(Tau_k,1) ~= 3 || size(Tau_k,2) ~= 3
            continue; % skip malformed facet
        end
        if size(U_k,1) ~= 3 || size(U_k,2) ~= N
            continue; % skip malformed facet
        end

        % Vertices as column vectors
        t1 = Tau_k(1,:).';  % 3x1
        t2 = Tau_k(2,:).';
        t3 = Tau_k(3,:).';

        % Build 3x2 basis from edges (t2 - t1, t3 - t1)
        B = [t2 - t1, t3 - t1];   % 3x2
        rhs = ad - t1;            % 3x1

        % Solve for [b;c] in ad = t1 + b*(t2-t1) + c*(t3-t1)
        % Skip facet if near-degenerate
        if rank(B) < 2
            continue;
        end
        bc = B \ rhs;  % 2x1
        b = bc(1);
        c = bc(2);
        a = 1 - b - c;

        % Inside-triangle test with tolerance
        if a >= -tol && b >= -tol && c >= -tol && a <= 1+tol && b <= 1+tol && c <= 1+tol
            % Clamp tiny tolerance excursions
            a = max(0, min(1, a));
            b = max(0, min(1, b));
            c = max(0, min(1, c));

            % Interpolate control: ud = a*u1 + b*u2 + c*u3
            u1 = U_k(1,:).';
            u2 = U_k(2,:).';
            u3 = U_k(3,:).';
            ud = a*u1 + b*u2 + c*u3;

            found = 1;
            index = k;
            weights = [a; b; c];
            return;
        end
    end

    % If we reach here, no facet contained ad (ud already zeroed)
    % You can handle infeasibility upstream (e.g., schedule across facets or report error).
end
