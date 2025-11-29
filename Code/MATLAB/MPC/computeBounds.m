function uBounds = computeBounds(A, u_min, u_max, mask)
% COMPUTEWRENCHBOUNDS  Compute per-axis box bounds on [Fx; Fy; tau]
% from A, u_min/u_max and a healthy mask.
%
% Inputs:
%   A      : 3 x m control effectiveness matrix
%   u_min  : m x 1 min thruster levels
%   u_max  : m x 1 max thruster levels
%   mask   : 1 x m logical, true = healthy, false = failed
%
% Output:
%   uBounds: struct with fields
%            Fx_min, Fx_max, Fy_min, Fy_max, Tau_min, Tau_max

    u_min = u_min(:);
    u_max = u_max(:);

    if nargin >= 4 && ~isempty(mask)
        % Failed thrusters: no authority
        u_min(~mask) = 0;
        u_max(~mask) = 0;
    end

    % extract rows
    rowFx  = A(1, :).';
    rowFy  = A(2, :).';
    rowTau = A(3, :).';

    % helper for one row
    [Fx_min, Fx_max]     = oneRowBounds(rowFx,  u_min, u_max);
    [Fy_min, Fy_max]     = oneRowBounds(rowFy,  u_min, u_max);
    [Tau_min, Tau_max]   = oneRowBounds(rowTau, u_min, u_max);

    uBounds.Fx_min  = Fx_min;
    uBounds.Fx_max  = Fx_max;
    uBounds.Fy_min  = Fy_min;
    uBounds.Fy_max  = Fy_max;
    uBounds.Tau_min = Tau_min;
    uBounds.Tau_max = Tau_max;
end


function [amin, amax] = oneRowBounds(a, u_min, u_max)
% linear function a^T u over box [u_min, u_max]

    amin = 0;
    amax = 0;

    for j = 1:length(a)
        aj = a(j);
        if aj > 0
            amax = amax + aj * u_max(j);
            amin = amin + aj * u_min(j);   %always zero
        elseif aj < 0
            amax = amax + aj * u_min(j);
            amin = amin + aj * u_max(j);
        else
            % aj == 0, contributes nothing
        end
    end
end
