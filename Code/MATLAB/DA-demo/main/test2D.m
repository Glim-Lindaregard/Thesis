function  test2D(Null)


%% --- Data setup ---

A = [1 -1  0;
     0  1  1];

% Full 3D control cube vertices
Uverts = [0 0 0 0 1 1 1 1;   % u1
          0 0 1 1 0 0 1 1;   % u2
          0 1 0 1 0 1 0 1];  % u3



removeMask = findMask(A,Uverts);

if nargin == 1
    removeMask = zeros(1,8);
end

Uverts(:, logical(removeMask)) = [];


visualize2D(Uverts,A);


function insideIdx = findMask(A,Uverts)

    XY = A*Uverts; %2x8

    % XY : 2×N matrix of mapped AMS points (A*Uverts)
% Returns logical row vector the same size as N: 1 = interior point

    % transpose so rows are points
    P = XY.';
    % build convex hull
    [K,~] = convhull(P(:,1), P(:,2));
    hullVerts = P(K(1:end-1), :);  % remove duplicate last vertex

    % Preallocate
    insideIdx = false(1, size(P,1));

    % compute signed distances to each hull edge
    % for point p inside, all signed distances should be < 0 (or > 0) consistently
    nEdges = size(hullVerts,1);
    for i = 1:size(P,1)
        p = P(i,:);
        s = zeros(1,nEdges);
        for e = 1:nEdges
            p1 = hullVerts(e,:);
            p2 = hullVerts(mod(e,nEdges)+1,:);
            edge = p2 - p1;
            normal = [edge(2), -edge(1)];  % outward normal (CCW)
            s(e) = dot(normal, p - p1);
        end
        % Point is interior if all s have same sign and none are zero
        if all(s < -1e-12) || all(s > 1e-12)
            insideIdx(i) = true;
        end
    end

end


function visualize2D(Uverts,A)
% simpleTest(removeMask)
%
% Visualizes the 3D control-space polytope (or subset)
% and its mapped 2D Attainable Moment Set (AMS).
%
% INPUT:
%   removeMask : 1x8 logical or numeric vector
%                1 = remove that vertex, 0 = keep
%                Example: simpleTest([1,0,1,0,0,0,0,0])


%% --- Map to 2D output (AMS projection onto z=0 plane) ---
XY = A * Uverts;             % 2 x N
X = XY(1,:);  Y = XY(2,:);
T3 = [X; Y; zeros(1,size(Uverts,2))];

% 2D convex hull (AMS outline)
[XYu,~,~] = unique([X(:) Y(:)], 'rows');
if size(XYu,1) >= 3
    K2 = convhull(XYu(:,1), XYu(:,2));
else
    K2 = 1:size(XYu,1);
end

%% --- 3D convex hull of remaining vertices ---
if size(Uverts,2) >= 4
    K3 = convhull(Uverts(1,:)', Uverts(2,:)', Uverts(3,:)');
else
    K3 = boundary(Uverts(1,:)', Uverts(2,:)', Uverts(3,:)', 1);
end

%% --- Plot setup ---
figure('Color','k'); hold on; axis equal; grid on; view(3);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w');
xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
title('Control Space (3D) and AMS Projection','Color','w');

%% --- Plot shaded control-space geometry ---
patch('Vertices',Uverts','Faces',K3, ...
      'FaceColor',[0.3 0.6 1.0], ...
      'FaceAlpha',1, ...
      'EdgeColor',[0.8 0.8 1], ...
      'LineWidth',1.2);

% Vertices
scatter3(Uverts(1,:), Uverts(2,:), Uverts(3,:), ...
         45, 'w', 'filled', 'MarkerEdgeColor','k');

%% --- Plot mapped 2D AMS points ---
plot3(T3(1,:), T3(2,:), T3(3,:), 'o', ...
      'MarkerFaceColor',[0.3 0.8 1], ...
      'MarkerEdgeColor','none', 'MarkerSize',6);

% AMS hull outline (cyan)
plot3(XYu(K2,1), XYu(K2,2), zeros(size(K2)), 'c-', 'LineWidth',2);

%% --- Mapping arrows (always from 3D vertex to AMS point, with arrowheads) ---
for k = 1:size(Uverts,2)
    u = Uverts(:,k);   % start (control-space vertex)
    a = T3(:,k);       % end (mapped AMS point)
    v = a - u;         % vector from u → a

    % Use quiver3 for arrowhead
    quiver3(u(1), u(2), u(3), ...
            v(1), v(2), v(3), ...
            0, ...                    % no autoscaling
            'Color', 'r', ...
            'LineWidth', 1.3, ...
            'MaxHeadSize', 0.4, ...   % arrowhead relative size
            'AutoScale', 'off');      % prevent scaling distortions
end


%% --- Null-space line in green (legend-safe) ---
n = null(A);
if ~isempty(n)
    n = n / norm(n);
    umin = min(Uverts,[],2); umax = max(Uverts,[],2);
    center = (umin + umax)/2;  scale = max(umax - umin);

    t = [-1.5, 1.5];
    L = center + scale * n * t;

    % main line used for the legend entry
    plot3(L(1,:), L(2,:), L(3,:), ...
                   'g-', 'LineWidth', 2, ...
                   'DisplayName','Null-space line');

    % arrowhead (hidden from legend)
    v = diff(L,1,2);
    quiver3(L(1,1), L(2,1), L(3,1), v(1), v(2), v(3), 0, ...
            'Color','g', 'LineWidth',2, 'MaxHeadSize',0.5, ...
            'AutoScale','off', 'HandleVisibility','off');
end



%% --- Legend ---
legend({'Control-space faces','Vertices','Mapped points', ...
        'AMS hull (2D)','Mapping arrows','Null-space line'}, ...
       'TextColor','w','Location','northeastoutside');
end

end