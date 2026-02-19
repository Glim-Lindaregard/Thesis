
function visualizeAMS(U, Asys, mode, aProduced, ad, abc)
%VISUALIZEAMS  Subplot-safe AMS visualization (no figure/axes side effects)

    % ---------------- Defaults ----------------
    if nargin < 3 || isempty(mode),       mode       = "normal"; end
    if nargin < 4 || isempty(aProduced),  aProduced  = [0;0;0]; end
    if nargin < 5 || isempty(ad),         ad         = [0;0;0]; end
    if nargin < 6 || isempty(abc),        abc        = [0;0;0]; end

    mode = lower(string(mode));
    doExport = (mode == "export");

    % ---------------- Styles ----------------
    switch mode
        case "export"
            FaceColor  = [0.8 0.95 0.95];
            EdgeColor  = [0.3 0.3 0.3];
            FontSize   = 20;
            ViewAngles = [35 20];
            ShowBasis  = false;
        otherwise  % "normal"
            FaceColor  = 0.75*[0.8 0.95 0.95];
            EdgeColor  = [0 0 0];
            FontSize   = 10;
            ViewAngles = [45 25];
            ShowBasis  = false;
    end

    FaceAlpha = 1.0;
    LineWidth = 1.5;
    ShowProduced = false;
    ShowDesired  = false;

    % ---------------- Axes setup ----------------
    ax = gca;
    % Dummy handles for legend
    hSimplified = patch(ax, ...
        'Faces',[], 'Vertices',[], ...
        'FaceColor',[0.2 0.6 0.9], ...
        'EdgeColor',[0.3 0.3 0.3], ...
        'FaceAlpha',0.8);
    
    hTrue = patch(ax, ...
        'Faces',[], 'Vertices',[], ...
        'FaceColor',[0.85 0.2 0.2], ...
        'EdgeColor',[0.3 0.3 0.3], ...
        'FaceAlpha',0.8);

    hold(ax,'on');
    %axis(ax,'equal');
    axis(ax,'vis3d');
    %view(ax, ViewAngles);
    %view(3);

    set(gcf,'Color','w');
    set(ax,'Color','w');
    %grid(ax,"minor")
    
    ax.GridColor      = 0.3*[1 1 1];
    ax.GridLineStyle = "--";
    ax.MinorGridColor = EdgeColor;
    ax.GridAlpha      = 0.8;
    ax.MinorGridAlpha = 0.8;
    
    ax.XColor = EdgeColor;
    ax.YColor = EdgeColor;
    ax.ZColor = EdgeColor;
    
    ax.LineWidth = 1.2;
    %ax.FontName = 'Times New Roman';

    scale = 2;
    xlim(ax,[-scale scale]);
    ylim(ax,[-scale scale]);
    zlim(ax,[-scale/4 scale/4]);

    % ---------------- Draw AMS facets ----------------
    nFacets = size(U,3);
    for k = 1:nFacets
        Vk = Asys * U(:,:,k);

        verts = Vk.';
        faces = [1 2 3 4];

        patch(ax, ...
            'Faces',faces, ...
            'Vertices',verts, ...
            'FaceColor',FaceColor, ...
            'FaceAlpha',FaceAlpha, ...
            'EdgeColor',EdgeColor, ...
            'LineWidth',LineWidth, ...
            'HandleVisibility','off');
    end


    % ---------------- Produced / desired ----------------
    origin = [0 0 0];

    if ShowProduced
        quiver3(ax, origin(1),origin(2),origin(3), ...
            aProduced(1),aProduced(2),aProduced(3), ...
            'Color','b','LineWidth',1.2,'AutoScale','off');
    end

    if ShowDesired
        quiver3(ax, origin(1),origin(2),origin(3), ...
            ad(1),ad(2),ad(3), ...
            'Color','r','LineWidth',1.5,'AutoScale','off');
    end

    % ---------------- Optional basis vectors ----------------
    if ShowBasis && ~isempty(U)
        Vi = Asys * U(:,:,1);
        adi = Vi(:,1);
        adj = Vi(:,2);
        adk = Vi(:,4);

        quiver3(ax,0,0,0,adi(1),adi(2),adi(3),'Color',[0 0.6 0],'LineWidth',1.2);
        quiver3(ax,0,0,0,adj(1),adj(2),adj(3),'Color',[0 0.3 0.7],'LineWidth',1.2);
        quiver3(ax,0,0,0,adk(1),adk(2),adk(3),'Color',[0.7 0.3 0],'LineWidth',1.2);
    end

        
    ShowConstraint = true;

    if ShowConstraint && ~isempty(abc) && all(isfinite(abc(:))) && all(abs(abc(:)) > 1e-12)
        Fx_cap  = abs(abc(1));
        Fy_cap  = abs(abc(2));
        Tau_cap = abs(abc(3));

        % resolution of the surface
        nE = 40;

        % make a unit sphere then scale it into an ellipsoid
        [Xs, Ys, Zs] = sphere(nE);
        Xe = Fx_cap  * Xs;
        Ye = Fy_cap  * Ys;
        Ze = Tau_cap * Zs;

        % draw ellipsoid surface
        surf(ax, Xe, Ye, Ze, ...
            'FaceColor','r', ...
            'EdgeColor','r', ...
            'FaceAlpha',0.12, ...
            'EdgeAlpha',0.12, ...
            'LineWidth',0.5, ...
            'HandleVisibility','off');

    end
    % ---------------- Labels ----------------
    ax.FontSize = FontSize;
    xlabel(ax,'F_x (N)');
    ylabel(ax,'F_y (N)');
    zlabel(ax,'\tau_z (N\cdot m)');
    title('Attainable Moment Set (AMS)', 'color', EdgeColor);

    lgd = legend(ax, ...
        [hSimplified, hTrue], ...
        {'Simplified AMS bounds', 'True AMS'}, ...
        'Location','northoutside', ...
        'Orientation','horizontal');

    lgd.Box = 'off';
    lgd.TextColor = 0.3*[1 1 1];
    lgd.FontSize = 15;
    lgd.ItemTokenSize = [30 20]; 
    lgd.Location = "southoutside";
    % Manual positioning (normalized figure units)
    lgd.Units = 'normalized';
    pos = lgd.Position;
    pos(2) = pos(2) + 0.075;   % move DOWN (reduce gap)
    lgd.Position = pos;


    legend off
    grid(ax,'on');
    box(ax,'on');

    % ---------------- Export only if asked ----------------
    if doExport
        exportgraphics(ax,'AMS.pdf','ContentType','vector');
    end
end
