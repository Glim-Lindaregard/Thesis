function AMS = buildAMS_null(A, umin, umax)
% BUILDAMS_NULL  AMS via null-space/partition (per paper's post row-space description).
% -------------------------------------------------------------------------
% INPUTS
%   A     : 3×m, rank 3 (control effectiveness)
%   umin  : m×1 lower bounds
%   umax  : m×1 upper bounds
%
% OUTPUT (minimal, row-space-free)
%   AMS.A, AMS.umin, AMS.umax
%   AMS.S0, AMS.F0, AMS.Xi  : a null-space basis from a well-conditioned 3×3
%   AMS.pairs               : all unordered actuator pairs (i,j)
%   AMS.facets(k) with fields:
%       .ij         : [i j]
%       .locktype   : 'max' or 'min'
%       .vertsU     : 4×m control vertices [A;B;C;D]
%       .Tau   : 4×3 moment vertices  [A;B;C;D]

[rows,m] = size(A);
assert(rows==3,'A must be 3×m');  
umin = umin(:);  umax = umax(:);
assert(isvector(umin)&&isvector(umax)&&numel(umin)==m&&numel(umax)==m);

% 1) Choose a well-conditioned 3×3 (for Xi only; not used to build facets)
S0 = bestCondTriple(A, 1:m);
F0 = setdiff(1:m,S0);

% 2) Build a null-space basis Xi (columns span ker(A))
Xi = nullspaceBasisFromTriple(A, S0);

% 3) Facet construction in null-space/partition form (no normals stored)
pairs = nchoosek(1:m,2);
K = size(pairs,1);
facets(2*K,1) = struct('ij',[],'locktype','','vertsU',[],'Tau',[]);
fidx = 0;

for p = 1:K
    i = pairs(p,1); 
    j = pairs(p,2);

    % pick an invertible basic triple S excluding i,j; choose best-conditioned
    S  = bestCondTripleExcluding(A, [i j]);
    % locking rules using det signs
    s  = detSigns(A, i, j);
    [base_max, base_min] = lockingBases(s, umin, umax, i, j);

    % f_max,ij
    [U4,T4] = build_vertices(base_max, i, j, umin, umax, A);
    fidx = fidx+1;
    facets(fidx).ij = [i j];
    facets(fidx).locktype = 'max';
    facets(fidx).vertsU = U4; 
    facets(fidx).Tau = T4;

    % f_min,ij
    [U4,T4] = build_vertices(base_min, i, j, umin, umax, A);
    fidx = fidx+1;
    facets(fidx).ij = [i j];
    facets(fidx).locktype = 'min';
    facets(fidx).vertsU = U4; 
    facets(fidx).Tau = T4;
end

AMS = struct('A',A,'umin',umin,'umax',umax, ...
             'S0',S0,'F0',F0,'Xi',Xi, ...
             'pairs',pairs,'facets',facets);
end

function Sbest = bestCondTriple(A, idx)
% Return the 3-index subset with minimal condition number from idx
if nargin < 2, idx = 1:size(A,2); end
triples = nchoosek(idx,3);
bestCond = Inf; Sbest = triples(1,:);
for r = 1:size(triples,1)
    S = triples(r,:);
    cnd = cond(A(:,S));
    if isfinite(cnd) && cnd < bestCond
        bestCond = cnd; 
        Sbest = S;
    end
end
end

function S = bestCondTripleExcluding(A, exclude)
% Choose best-conditioned invertible 3×3 not using indices in 'exclude'
m = size(A,2);
cand = setdiff(1:m, exclude);
S = bestCondTriple(A, cand);
end

function Xi = nullspaceBasisFromTriple(A, S0)
% Build a null-space basis using a fixed well-conditioned S0
m = size(A,2);
F0 = setdiff(1:m, S0);
AS0 = A(:,S0);
Xi = zeros(m, numel(F0));
for k = 1:numel(F0)
    f = F0(k);
    xiS = -(AS0 \ A(:,f));  % solve AS0*xiS + A(:,f) = 0
    v = zeros(m,1); 
    v(S0) = xiS; 
    v(f)  = 1;
    Xi(:,k) = v;
end
end

function s = detSigns(A, i, j, tol)
% s(k) = sign( det([c_i c_j c_k]) ), with s(i)=s(j)=0
if nargin < 4, tol = 1e-12; end
m = size(A,2);
ci = A(:,i); 
cj = A(:,j);
s = zeros(m,1);
for k = 1:m
    if k==i || k==j
        s(k) = 0;
    else
        M = [ci cj A(:,k)];
        d = det(M);
        if abs(d) < tol, d = 0; end
        s(k) = sign(d);
    end
end
end

function [base_max, base_min] = lockingBases(s, umin, umax, i, j)
% Apply Ω_max,ij and Ω_min,ij rules; i,j remain free (not set here)
base_max = umin; 
base_max(s > 0) = umax(s > 0);

base_min = umin; 
base_min(s < 0) = umax(s < 0);

% ensure i,j are not locked by accident (leave as NaN to emphasize "free")
base_max([i j]) = NaN;
base_min([i j]) = NaN;
end


function [U4, T4] = build_vertices(base, i, j, umin, umax, A)
% Create the 4 vertices for free pair (i,j) with others fixed by 'base'
% U4: 4×m, rows are vertices; T4: 4×3, corresponding moments
m = numel(base);
choices = [umin(i) umin(j);
           umin(i) umax(j);
           umax(i) umin(j);
           umax(i) umax(j)];    % [ (i,j) = (min,min); (min,max); (max,min); (max,max) ]

U4 = zeros(4,m);
for r = 1:4
    u = base;
    % fill free indices
    u(i) = choices(r,1);
    u(j) = choices(r,2);
    % if any leftover NaNs (shouldn't happen), default to umin
    nanidx = isnan(u);
    if any(nanidx), u(nanidx) = umin(nanidx); end
    U4(r,:) = u.';
end

T4 = (A * U4.').';  % 4×3
end
