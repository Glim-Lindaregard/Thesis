S = load('planes_mode0.mat');
n = S.n;            % (P x 3)
b = S.b(:);         % (P x 1)
P = size(n, 1);

% Plot limits (adjust to taste)
Fx_lim_g  = [-2  2];
Fy_lim_g  = [-2  2];
Fx_lim = [-1.5, 1.5];
Fy_lim = [-1.5, 1.5];
Tau_lim_g = [-0.5  0.5];
Tau_lim = [-0.4, 0.4];


buff = 0.2;

% Grid resolution
L = 1.7;          % your quadrant extent
eps0 = 1e-8;      % threshold for "component is basically zero"
Ng = 50;

figure; hold on; grid on;
%xlabel('F_x'); ylabel('F_y'); zlabel('\tau');

%view(3);

% Plot origin
plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'k');

for i = 1:P
    ni = n(i, :);
    bi = b(i);

    % if ni(1) >= 0 && ni(2) >= 0
    %     Fx_lim = [0 Fl];
    %     Fy_lim = [0 Fl];
    % elseif ni(1) >= 0 && ni(2) <= 0 
    %     Fx_lim = [0 Fl];
    %     Fy_lim = [-Fl 0];
    % elseif ni(1) <= 0 && ni(2) <= 0
    %     Fx_lim = [-Fl 0]; 
    %     Fy_lim = [-Fl 0];
    % elseif ni(1) <= 0 && ni(2) >= 0
    %     Fx_lim = [-Fl 0]; 
    %     Fy_lim = [0 Fl];
    % end

    % Choose which variable to solve for
    [~, j] = max(abs(ni));  % 1->Fx, 2->Fy, 3->Tau

    switch j
        case 3  % solve for Tau
            if abs(ni(3)) < 1e-10, continue;  end
            [Fx_lim_q, Fy_lim_q] = quadrant_limits(ni(1), ni(2), L, eps0);
            [Fx, Fy] = meshgrid( ...
                linspace(Fx_lim_q(1), Fx_lim_q(2), Ng), ...
                linspace(Fy_lim_q(1), Fy_lim_q(2), Ng) ...
            );
            
            Tau = (bi - ni(1)*Fx - ni(2)*Fy) / ni(3) + [0 0 0.01]*ni';
            
            surf(Fx, Fy, Tau, ...
                'FaceColor', [0.2 0.6 0.9], ...
                'FaceAlpha', 0.5, ...
                'EdgeAlpha', 0.7, ...
                'EdgeColor', 0.5*[1 1 1]);

        case 2  % solve for Fy
            if abs(ni(2)) < 1e-10, continue;  end
            [Fx_lim_q, Fy_lim_q] = quadrant_limits(ni(1), ni(3), L, eps0);
            [Fx, Tau] = meshgrid( ...
                linspace(Fx_lim_q(1), Fx_lim_q(2), Ng), ...
                linspace(Fy_lim_q(1), Fy_lim_q(2), Ng) ...
            );
            Fy = (bi - ni(1)*Fx - ni(3)*Tau) / ni(2) + [0 0.01 0]*ni';

            surf(Fx, Fy, Tau, ...
                'FaceColor', [0.2 0.6 0.9], ...
                'FaceAlpha', 0.5, ...
                'EdgeAlpha', 0.7, ...
                'EdgeColor', 0.5*[1 1 1]);

        case 1  % solve for Fx
            if abs(ni(1)) < 1e-10, continue;  end
            [Fx_lim_q, Fy_lim_q] = quadrant_limits(ni(2), ni(2), L, eps0);
            [Fy, Tau] = meshgrid( ...
                linspace(Fx_lim_q(1), Fx_lim_q(2), Ng), ...
                linspace(Fy_lim_q(1), Fy_lim_q(2), Ng) ...
            );

            Fx = (bi - ni(2)*Fy - ni(3)*Tau) / ni(1) + [0.01 0 0]*ni';

            surf(Fx, Fy, Tau, ...
                'FaceColor', [0.2 0.6 0.9], ...
                'FaceAlpha', 0.5, ...
                'EdgeAlpha', 0.7, ...
                'EdgeColor', 0.5*[1 1 1]);
    end
end

c = config();
visualizeAMS(buildAMS_row(c), c.A, 'export');

axis([Fx_lim_g Fy_lim_g Tau_lim_g]);
%axis vis3d;
function [limA, limB] = quadrant_limits(na, nb, L, eps0)
% Returns limits for two axes based on signs of na, nb.
% If one component is near zero, use symmetric limits for that axis.
    if abs(na) < eps0
        limA = [-L L];
    elseif na >= 0
        limA = [0 L];
    else
        limA = [-L 0];
    end

    if abs(nb) < eps0
        limB = [-L L];
    elseif nb >= 0
        limB = [0 L];
    else
        limB = [-L 0];
    end
end
