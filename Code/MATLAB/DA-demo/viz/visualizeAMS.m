function visualizeAMS(AMS, opts)
% Plot the AMS facets with configurable "Tang-style" aesthetics.
%
% opts (all optional):
%   Appearance
%     .FaceColor        (default [0.78 0.92 0.92])  % pastel cyan
%     .EdgeColor        (default 'k')               % black edges
%     .FaceAlpha        (default 0.95)
%     .EdgeAlpha        (default 0.85)
%     .LineWidth        (default 0.75)
%     .UseOctantColors  (default false)             % override FaceColor per facet
%     .BackgroundColor  (default 'w')
%     .Lighting         (default false)             % add a soft camlight
%
%   Axes / grid
%     .GridColor        (default [0.55 0.55 0.55])
%     .GridAlpha        (default 0.35)
%     .GridLineStyle    (default '--')
%     .FontSize         (default 11)
%     .AxisLabels       (default {'T_{ex} (N·m)','T_{ey} (N·m)','T_{ez} (N·m)'})
%     .AxisLimits       (default [])  % e.g., [-1 1; -1 1; -0.5 0.5]
%   Geometry
%     .ShowNormals      (default false)
%     .MaxFacets        (default inf)


% --- Options ---
    if nargin<2, opts = struct(); end

    FaceColor       = getfield_with_default(opts,'FaceColor',[0.78 0.92 0.92]);
    EdgeColor       = getfield_with_default(opts,'EdgeColor','k');
    FaceAlpha       = getfield_with_default(opts,'FaceAlpha',0.95);
    EdgeAlpha       = getfield_with_default(opts,'EdgeAlpha',0.85);
    LineWidth       = getfield_with_default(opts,'LineWidth',0.75);
    BackgroundColor = getfield_with_default(opts,'BackgroundColor','w');
    LightingOn      = getfield_with_default(opts,'Lighting',false);

    GridColor     = getfield_with_default(opts,'GridColor',[0.55 0.55 0.55]);
    GridAlpha     = getfield_with_default(opts,'GridAlpha',0.35);
    GridStyle     = getfield_with_default(opts,'GridLineStyle','--');
    FontSize      = getfield_with_default(opts,'FontSize',11);
    AxisLimits    = getfield_with_default(opts,'AxisLimits',[]);

    ShowNormals   = getfield_with_default(opts,'ShowNormals',false);
    MaxFacets     = getfield_with_default(opts,'MaxFacets',inf);



    % --- Figure --- 

    fig = figure('Color',BackgroundColor); %#ok<NASGU>
    hold on;

    count = min(numel(AMS.facets), MaxFacets);

    for k = 1:count
        Tau = AMS.facets(k).Tau;   % 4x3, order A,B,C,D

        colFace = FaceColor;
        colEdge = EdgeColor;

        % --- Split facets into two triangles for better rendering --- 
        A = Tau(1,:); B = Tau(2,:); C = Tau(3,:); D = Tau(4,:);
        tris = [A; B; C; A; C; D];
        f = [1 2 3; 4 5 6];

        patch('Faces',f,'Vertices',tris, ...
              'FaceColor',colFace, 'FaceAlpha',FaceAlpha, ...
              'EdgeColor','none', 'EdgeAlpha',EdgeAlpha, ...
              'LineWidth',LineWidth, ...
              'EdgeLighting','none','FaceLighting','flat');

        outline = [A;B;C;D;A];

        plot3(outline(:,1), outline(:,2), outline(:,3), ...
              'Color', colEdge, 'LineWidth', LineWidth, ...
              'LineJoin','round');
        
        % --- generate normals --- 
        if ShowNormals
            ctr = mean(Tau,1);
            n   = AMS.facets(k).normal(:);
            if norm(n) > 0
                n = n/norm(n);
                quiver3(ctr(1),ctr(2),ctr(3), 0.3*n(1),0.3*n(2),0.3*n(3), ...
                        'AutoScale','off','Color','r','LineWidth',2, 'MaxHeadSize', 1);
            end
        end
    end

    % Axes styling
    axis equal; box on; grid on;
    ax = gca;
    ax.GridColor     = GridColor;
    ax.GridAlpha     = GridAlpha;
    ax.GridLineStyle = GridStyle;
    ax.MinorGridLineStyle = GridStyle;
    ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on'; ax.ZMinorGrid = 'on';
    ax.FontSize = FontSize;

    % --- Supposed to be labes editable from opts but it does not work (switches to default)---
    lbl = getfield_with_default(opts,'AxisLabels', ...
       {'F_{x} (N)','F_{y} (N)','T_{z} (N·m)'});

    xlabel(lbl(1));
    ylabel(lbl(2));
    zlabel(lbl(3));
    title('Attainable Moment Set');



    if ~isempty(AxisLimits)
        ax.XLim = AxisLimits(1,:);
        ax.YLim = AxisLimits(2,:);
        ax.ZLim = AxisLimits(3,:);
    end


    view(35,25);

    % --- Dotted origin axes (subtle) ---
    lim = [ax.XLim; ax.YLim; ax.ZLim];
    plot3([0 0], lim(2,:), [0 0], ':', 'Color', [0 0 0], 'LineWidth', 0.5);
    plot3(lim(1,:), [0 0], [0 0], ':', 'Color', [0 0 0], 'LineWidth', 0.5);
    plot3([0 0], [0 0], lim(3,:), ':', 'Color', [0 0 0], 'LineWidth', 0.5);

    % --- Optional soft lighting ---
    if LightingOn
        camlight('headlight'); material dull; lighting gouraud; %#ok<UNRCH>
    end

end

% function idx = octant_index(oct)
%     if any(oct==0), idx = 9; return; end
%     b = (oct>0);
%     map = [1;2;3;4;5;6;7;8];
%     idx = map( bin2dec(sprintf('%d%d%d', b(1), b(2), b(3))) + 1 );
% end


% --- Function used to input options ---
function v = getfield_with_default(s, name, def)
    if isfield(s,name), v = s.(name); else, v = def; end
end
