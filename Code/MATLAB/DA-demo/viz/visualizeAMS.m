function visualizeAMS(U, Asys, mode, aProduced, ad, abc)

    % --- defaults for optional vectors ---
    if nargin < 4 || isempty(aProduced), aProduced = [0;0;0]; end
    if nargin < 5 || isempty(ad),        ad        = [0;0;0]; end
    if nargin < 6 || isempty(abc),      abc       = [0;0;0]; end


    mode = lower(string(mode));
    doExport = strcmp(mode, "export");

    % --- preset styles ---
    switch mode
        case "export"
            FaceColor       = [0.8 0.95 0.95];
            EdgeColor       = [0.3 0.3 0.3];
            FaceAlpha       = 0.95;
            EdgeAlpha       = 1;
            LineWidth       = 0.8;
            BackgroundColor = [1 1 1];
            GridColor       = 0.3*[1 1 1];
            GridAlpha       = 1;
            GridStyle       = '--';
            FontSize        = 30;
            ViewAngles      = [55 30];
            ShowProduced    = true;
            ShowDesired     = true;
            ShowBasis       = true;
        case "normal"
            FaceColor       = 0.75*[0.8 0.95 0.95];
            EdgeColor       = [0 0 0];
            FaceAlpha       = 0.95;
            EdgeAlpha       = 1;
            LineWidth       = 0.8;
            BackgroundColor = 0.2*[1 1 1];
            GridColor       = 0.7*[1 1 1];
            GridAlpha       = 1;
            GridStyle       = '--';
            FontSize        = 11;
            ViewAngles      = [55 25];
            ShowProduced    = true;
            ShowDesired     = true;
            ShowBasis       = false;
        otherwise
            error('visualizeAMS:UnknownMode','Mode must be "normal" or "export".');
    end

    facetIndex = 1;   % which facet to highlight

    % --- basic figure setup ---
    figure('Color', BackgroundColor); hold on;
    set(gcf,'Renderer','painters');

    count  = size(U,3);
    scale  = 2;
    center = [0 0 0];

    axis manual; axis equal; axis vis3d;
    xlim([-scale,scale]);
    ylim([-scale,scale]);
    zlim([-scale/4,scale/4]);
    view(ViewAngles(1), ViewAngles(2));

    lgHandles = [];
    lgLabels  = {};

    % --- draw facets ---
    for k = 1:count
        Uk   = U(:,:,k);
        Vk   = Asys*Uk;
        verts = Vk;

        A = verts(:,1)'; B = verts(:,2)'; C = verts(:,3)'; D = verts(:,4)';
        tris = [A; B; C; D];
        f    = [1 2 3 4];

        isHighlighted = (k == facetIndex) && (ShowBasis || ShowProduced);

        if isHighlighted
            patch('Faces',f,'Vertices',tris, ...
              'FaceColor','r','FaceAlpha',0.5, ...
              'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
              'EdgeLighting','none','FaceLighting','flat','HandleVisibility','off');
        else
            patch('Faces',f,'Vertices',tris, ...
              'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
              'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
              'EdgeLighting','none','FaceLighting','flat','HandleVisibility','off');
        end

        % highlight edges of selected facet
        if isHighlighted
            outline = [A;B;C;D;A];
            for i = 1:4
                plot3(outline(i:i+1,1),outline(i:i+1,2),outline(i:i+1,3), ...
                      'Color',[0,i*0.24,0], ...
                      'LineWidth',LineWidth,'LineJoin','round');
            end
        end

    end

    % --- produced / desired wrenches ---
    if ShowProduced
        hProd = quiver3(center(1),center(2),center(3), ...
                        aProduced(1), aProduced(2), aProduced(3), ...
                        'AutoScale','off','Color','b', ...
                        'LineWidth',1.5,'MaxHeadSize',0.8);
        lgHandles(end+1) = hProd;
        lgLabels{end+1}  = 'Moment produced';
    end

    if ShowDesired
        hDes = quiver3(center(1),center(2),center(3), ...
                       ad(1), ad(2), ad(3), ...
                       'AutoScale','off','Color','r', ...
                       'LineWidth',2,'MaxHeadSize',0.8);
        lgHandles(end+1) = hDes;
        lgLabels{end+1}  = 'Desired moment';
    end

    % --- basis vectors for selected facet ---
    if ShowBasis
        Vi  = Asys*U(:,:,facetIndex);
        adi = Vi(:,1);
        adj = Vi(:,2);
        adk = Vi(:,4);

        a = abc(1); b = abc(2); c = abc(3);

        M  = [adi, b*(adj-adi), c*(adk-adi)];
        M2 = [adi-center' adj adk];

        for n = 1:3
            h1 = quiver3(center(1),center(2),center(3), ...
                         M2(1,n), M2(2,n), M2(3,n), ...
                         'AutoScale','off','Color',[0,n*0.33,n*0.33], ...
                         'LineWidth',1.5,'MaxHeadSize',0.8);
            h2 = quiver3(adi(1),adi(2),adi(3), ...
                         M(1,n), M(2,n), M(3,n), ...
                         'AutoScale','off','Color',[1,0.22,0], ...
                         'LineWidth',1.5,'MaxHeadSize',0.8);

            if n == 3
                lgHandles(end+1) = h1;
                lgLabels{end+1}  = 'adi adj adk';
                lgHandles(end+1) = h2;
                lgLabels{end+1}  = 'Facet basis vectors';
            end
        end
    end

    % --- axes / labels / legend ---
    box on; grid on;
    ax = gca;
    ax.GridColor     = GridColor;
    ax.GridAlpha     = GridAlpha;
    ax.GridLineStyle = GridStyle;
    ax.XMinorGrid    = 'off';
    ax.YMinorGrid    = 'off';
    ax.ZMinorGrid    = 'off';
    ax.FontSize      = FontSize;
    ax.LineWidth     = LineWidth;
    ax.Color         = 0.1*[1 1 1];
    ax.Projection    = 'perspective';
    ax.XColor        = GridColor;
    ax.YColor        = GridColor;
    ax.ZColor        = GridColor;

    xlabel('F_x (N)','Interpreter','tex');
    ylabel('F_y (N)','Interpreter','tex');
    zlabel('T_z (N·m)','Interpreter','tex');
    title('Attainable Moment Set','Color',GridColor);

    if ~isempty(lgHandles)
        legend(lgHandles, lgLabels{:});
    end

    % only export pdf in export mode
    if doExport
        exportgraphics(gcf, 'AMS.pdf', 'ContentType','vector','Padding',5);
    end
end
