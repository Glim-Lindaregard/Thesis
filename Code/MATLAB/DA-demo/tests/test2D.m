function  test2D()
close all;

% --- Data setup ---
A = [1 -1  0;
     0  -1  1];

% Full 3D control cube vertices
Uverts = [0 0 0 0 1 1 1 1;   % u1
          0 0 1 1 0 0 1 1;   % u2
          0 1 0 1 0 1 0 1];  % u3

ad = [-2 -1]';

%visualize2D(Uverts,A);
search2D(ad,Uverts,A);


function visualize2D(Uverts,A)

%%% --- Map to 2D output (AMS projection onto z=0 plane) ---
XY = A * Uverts;             % 2 x N
X = XY(1,:);  Y = XY(2,:);
T3 = [X; Y; zeros(1,size(Uverts,2))];

%2D convex hull (AMS outline)
K2 = convhull(X, Y);

% --- 3D convex hull of remaining vertices ---
K3 = convhull(Uverts(1,:)', Uverts(2,:)', Uverts(3,:)');

% --- Plot setup ---
figure('Color','k'); hold on; axis equal; grid on; view(3);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w');
xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
title('Control Space (3D) and AMS Projection','Color','w');

% --- Plot shaded control-space geometry ---
trisurf(K3, Uverts(1,:), Uverts(2,:), Uverts(3,:), ...
    'FaceColor', [0.3 0.6 1.0], ...   % solid blueish color
    'EdgeColor', [0.8 0.8 1], ...     % soft white edges
    'LineWidth', 1.2, ...
    'FaceAlpha',0.5);

% Vertices
scatter3(Uverts(1,:), Uverts(2,:), Uverts(3,:), ...
         45, 'w', 'filled', 'MarkerEdgeColor','k');

% --- Plot mapped 2D AMS points ---
plot3(T3(1,:), T3(2,:), T3(3,:), 'o', ...
      'MarkerFaceColor',[0.3 0.8 1], ...
      'MarkerEdgeColor','none', 'MarkerSize',6);

%AMS hull outline (cyan)
plot3(X(K2), Y(K2), zeros(size(K2)), 'c-', 'LineWidth',2);
% --- Mapping arrows (always from 3D vertex to AMS point, with arrowheads) ---
for k = 1:size(Uverts,2)
    u = Uverts(:,k);   % start (control-space vertex)
    a = T3(:,k);       % end (mapped AMS point)
    v = a - u;         % vector from u → a

    % Use quiver3 for arrowhead
    quiver3(u(1), u(2), u(3), ...
            v(1), v(2), v(3), ...
            0, ...                    % no autoscaling
            'Color', [1 0.2 0.2], ...
            'LineWidth', 1.5, ...
            'MaxHeadSize', 0.4, ...   % arrowhead relative size
            'AutoScale', 'off');      % prevent scaling distortions
end





% --- Legend ---
legend({'Control-space faces','Vertices', ...
       'AMS hull (2D)'}, ...
      'TextColor','w','Location','northeastoutside');

end

function search2D(ad,Uverts,A)
%%% --- Map to 2D output (AMS projection onto z=0 plane) ---
XY = A * Uverts;             % 2 x N
X = XY(1,:);  Y = XY(2,:);

%2D convex hull (AMS outline)
K2 = convhull(X, Y);

bestLam = inf; bestPt = [NaN; NaN];
for k = 1:length(K2)-1
    ai = XY(:, K2(k));
    aj = XY(:, K2(k+1));
    M  = [ad, -(aj - ai)];
    if abs(det(M)) < 1e-9, continue; end
    x = M \ ai;                % [lambda; alpha]
    lam = x(1); alpha = x(2);
    if lam > 0 && alpha >= 0 && alpha <= 1 && lam < bestLam
        bestLam = lam;
        bestPt  = lam * ad;   % intersection point on AMS
    end
end

% --- Plot setup ---
figure('Color','k'); hold on; axis equal; grid on; view(3);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w');
xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
title('Control Space (3D) and AMS Projection','Color','w');


%AMS hull outline (cyan)
plot3(X(K2), Y(K2), zeros(size(K2)), 'c-', 'LineWidth',2);

plot3([0 ad(1)], [0 ad(2)], [0 0], 'm--', 'LineWidth', 1.5);
plot3([0 bestPt(1)], [0 bestPt(2)], [0 0], 'm-', 'LineWidth', 2);
scatter3(bestPt(1), bestPt(2), 0, 80, 'm', 'filled');

legend('AMS','desired','produced','intersection');

end


end
