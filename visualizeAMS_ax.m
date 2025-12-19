function visualizeAMS_ax(ax, U, Asys, mode)
% Plot AMS facets into a given axes handle (no figure creation, no export).
% Expects U as N x 4 x F.

mode = lower(string(mode));

% --- preset styles (export tuned) ---
switch mode
    case "export"
        FaceColor   = [0.8 0.95 0.95];
        EdgeColor   = 0.5*[0.3 0.3 0.3];
        FaceAlpha   = 0.95;
        EdgeAlpha   = 1;
        LineWidth   = 0.8;
        GridColor   = 0.2*[1 1 1];
        GridAlpha   = 1;
        GridStyle   = '-';
        FontSize    = 10;      % smaller for 2x2 page
        ViewAngles  = [55 30];
        scale       = 2.6;     % zoom OUT
    otherwise
        FaceColor   = 0.75*[0.8 0.95 0.95];
        EdgeColor   = [0 0 0];
        FaceAlpha   = 0.95;
        EdgeAlpha   = 1;
        LineWidth   = 0.8;
        GridColor   = 0.7*[1 1 1];
        GridAlpha   = 1;
        GridStyle   = '--';
        FontSize    = 11;
        ViewAngles  = [55 25];
        scale       = 2.3;
end

% --- axes setup ---
cla(ax);
hold(ax,'on');
axis(ax,'equal'); axis(ax,'vis3d');
set(ax,'Projection','perspective');

xlim(ax, [-scale,  scale]);
ylim(ax, [-scale,  scale]);
zlim(ax, [-scale/4, scale/4]);
view(ax, ViewAngles(1), ViewAngles(2));

% WHITE axis walls, white legend box, black text
ax.Color          = [1 1 1];     % axis walls
ax.XColor         = GridColor;
ax.YColor         = GridColor;
ax.ZColor         = GridColor;
ax.GridColor      = GridColor;
ax.GridAlpha      = GridAlpha;
ax.GridLineStyle  = GridStyle;
ax.FontSize       = FontSize;
ax.LineWidth      = LineWidth;

grid(ax,'on'); box(ax,'on');

xlabel(ax,'F_x (N)','Interpreter','tex');
ylabel(ax,'F_y (N)','Interpreter','tex');
zlabel(ax,'T_z (N·m)','Interpreter','tex');

% --- draw facets ---
count = size(U,3);
for k = 1:count
    Uk = U(:,:,k);
    Vk = Asys * Uk;   % 3x4

    A1 = Vk(:,1)'; B1 = Vk(:,2)'; C1 = Vk(:,3)'; D1 = Vk(:,4)';
    tris = [A1; B1; C1; D1];
    f = [1 2 3 4];

    patch(ax,'Faces',f,'Vertices',tris, ...
        'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
        'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
        'EdgeLighting','none','FaceLighting','flat','HandleVisibility','off');
end




% no basis/desired/produced => no legend should exist, but force style anyway
% lg = legend('axis');
% if isvalid(lg)
%     lg.Color     = [1 1 1];
%     lg.TextColor = [0 0 0];
% end
end
