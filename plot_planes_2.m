S = load('planes_mode0.mat');
n = S.n;  b = S.b(:);

% ---------- Plot ranges / resolution ----------
Fx = linspace(-1.6, 1.6, 100);
Fy = linspace(-1.6, 1.6, 100);
Tau = linspace(-0.5, 0.5, 100);

tol = 1e-7;        % inequality tolerance for clipping
nz_tol = 1e-6;     % "vertical plane" threshold
% --------------------------------------------

cfg0 = config();
uCache = buildAllAMS(cfg0);  



figure; hold on; grid on; axis equal
xlabel('Fx'); ylabel('Fy'); zlabel('Tau');
title('AMS planes clipped to polytope (grid faces only)');
plot3(0,0,0,'k.','MarkerSize',18);

% helper: check all halfspaces n*x <= b
% X,Y,Z are same size arrays
feasible_mask = @(X,Y,Z) all( ...
    reshape( n(:,1).*X(:).' + n(:,2).*Y(:).' + n(:,3).*Z(:).', size(n,1), [] ) ...
    <= (b + tol), ...
    1);

for i = 1:size(n,1)
    nx = n(i,1); ny = n(i,2); nz = n(i,3);

    col = 0.5*(n(i,:)+1); % color by normal direction

    if abs(nz) > nz_tol
        % ---- regular: Tau = (b - nx*Fx - ny*Fy)/nz ----
        [FX,FY] = meshgrid(Fx, Fy);
        TAU = (b(i) - nx*FX - ny*FY) / nz;

        mask = feasible_mask(FX, FY, TAU);
        mask = reshape(mask, size(FX));

        TAU(~mask) = NaN; % clip to polytope

        mesh(FX, FY, TAU, ...
            'EdgeColor', col, ...
            'FaceColor', 'none', ...
            'LineWidth', 1.0);

    else
        % ---- vertical-ish: nx*Fx + ny*Fy = b ----
        if abs(nx) >= abs(ny) && abs(nx) > 1e-9
            % Solve for Fx, sweep Fy and Tau
            [FY2,TAU2] = meshgrid(Fy, Tau);
            FX2 = (b(i) - ny*FY2) / nx;

            mask = feasible_mask(FX2, FY2, TAU2);
            mask = reshape(mask, size(FX2));

            FX2(~mask) = NaN;

            mesh(FX2, FY2, TAU2, ...
                'EdgeColor', col, ...
                'FaceColor', 'none', ...
                'LineWidth', 1.0);

        elseif abs(ny) > 1e-9
            % Solve for Fy, sweep Fx and Tau
            [FX2,TAU2] = meshgrid(Fx, Tau);
            FY2 = (b(i) - nx*FX2) / ny;

            mask = feasible_mask(FX2, FY2, TAU2);
            mask = reshape(mask, size(FY2));

            FY2(~mask) = NaN;

            mesh(FX2, FY2, TAU2, ...
                'EdgeColor', col, ...
                'FaceColor', 'none', ...
                'LineWidth', 1.0);
        end
        hold on;
        health_mask = [1 1 1 1 1 1 1 1]; %%%%%%%%%%%%%%%%%%%%%%
        idx = mask_to_index(health_mask);
        if any(health_mask==2)
            idx2 = find(health_mask == 2);
            idx0 = find(health_mask == 0);
            cfg0.u_min(idx2) = 0.7;
            cfg0.u_max(idx0) = 0.0;
            U = buildAMS_row(cfg0);
            hold on;
            visualizeAMS(U, cfg0.A, 'normal');
        else
            visualizeAMS(uCache(idx+1).U,uCache(idx+1).A,'normal');
        end
    end
end

view(35,25)

function idx = mask_to_index(mask)
    idx = 0;
    for k = 1:numel(mask)
        if mask(k) ~= 0
            idx = bitor(idx, bitshift(1, k-1));
        end
    end
end
function uCache = buildAllAMS(cfg0)
    % --- Base config ---
    N     = cfg0.N_thrusters;                 % should be 8
     
    %Cashe to store all 256 AMSs and As
    uCache = struct([]);
    
    %Build all 256 AMSs
    for k = 0:(2^N - 1)
        mask = logical(bitget(k, 1:N));   % 1 = healthy, 0 = failed
        
        cfgi = cfg0;
        cfgi.A(:, ~mask) = 0;             % zero failed thruster columns
        cfgi.u_min(~mask) = 0;            % zero min/max ranges
        cfgi.u_max(~mask) = 0;
        cfgi.N = N;                       % keep same dimension
        
        Ui= buildAMS_row(cfgi);
        uCache(k+1).U = Ui;
        uCache(k+1).A = cfgi.A;
        uCache(k+1).mask = mask;
    end
end