function visualizeAMS(AMS, opts,aProduced,ad,abc)

%   --- Options ---
if nargin<2, opts=struct(); end
FaceColor       = getfield_with_default(opts,'FaceColor',[0.78 0.92 0.92]);
EdgeColor       = getfield_with_default(opts,'EdgeColor','k');
FaceAlpha       = getfield_with_default(opts,'FaceAlpha',0.95);
EdgeAlpha       = getfield_with_default(opts,'EdgeAlpha',0.85);
LineWidth       = getfield_with_default(opts,'LineWidth',0.75);
BackgroundColor = getfield_with_default(opts,'BackgroundColor','w');
LightingOn      = getfield_with_default(opts,'Lighting',false);
GridColor       = getfield_with_default(opts,'GridColor',[0.55 0.55 0.55]);
GridAlpha       = getfield_with_default(opts,'GridAlpha',0.35);
GridStyle       = getfield_with_default(opts,'GridLineStyle','--');
FontSize        = getfield_with_default(opts,'FontSize',11);
ShowNormals     = getfield_with_default(opts,'ShowNormals',false);
fps             = getfield_with_default(opts,'fps',10);
fancy           = getfield_with_default(opts,'Fancy',true);
index           = getfield_with_default(opts,'Index',0);
ShowProduced        = getfield_with_default(opts,'ShowProduced',true);
ShowDesired         = getfield_with_default(opts,'ShowDesired',true);
ShowBasis           = getfield_with_default(opts,'ShowBasis',true);

figure('Color',BackgroundColor); hold on
count = numel(AMS.facets);

scale = 2;

axis manual; axis equal; axis vis3d; ...
    xlim([-scale,scale]); ylim([-scale,scale]); zlim([-scale/2,scale/2]);
view(45,30) 

center = [0 0 0]; %AMS.center;

for k = 1:count
    %Get a facets verteces.
    verts = AMS.facets(k).V;
    
    %extract verteces
    A = verts(:,1)'; B = verts(:,2)'; C = verts(:,3)'; D = verts(:,4)';
    
    %tris = [A; B; C; A; C; D];  If plotting looks weird, might be this.

    tris = [A; B; C; D];
    
    %f = [1 2 3; 4 5 6];   And this

    f = [1 2 3 4];

    %Plot facets
    if k == index
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

    %Fancy camera spin
    if fancy
        for i = 1:fps
            camorbit(360/(fps*count),0,'data'); 
            drawnow
            %pause(0.002)
        end
        drawnow limitrate
    end


    %Make selected facet edges red
    if k == index
        for i = 1:4
            %Edges %Put back on if edges look wierd
            outline = [A;B;C;D;A]; 
            plot3(outline(i:i+1,1),outline(i:i+1,2),outline(i:i+1,3), ...
                  'Color',[0,i*0.24,0],'LineWidth',LineWidth,'LineJoin','round');
        end
    end
    
    %Show normals
    if ShowNormals
        n = AMS.facets(k).norms;
        ctr = mean(verts,2);
        %Make sure normals point in "correct" direction
        if dot(ctr-center,n)<0
            n = -n;
        end
        quiver3(ctr(1),ctr(2),ctr(3),0.2*n(1),0.2*n(2),0.2*n(3), ...
        'AutoScale','off','Color','r','LineWidth',1.5,'MaxHeadSize',...
        0.8,'HandleVisibility','off');
    end
end

if ShowProduced
    producedHandle = quiver3(center(1),center(2),center(3),1*aProduced(1),1*aProduced(2),1*aProduced(3)...
    ,'off','Color','b','LineWidth',1.5,'MaxHeadSize',0.8);
end

if ShowDesired
    desiredHandle = quiver3(center(1),center(2),center(3),1*ad(1),1*ad(2),1*ad(3)...
    ,'off','Color','g','LineWidth',0.5,'MaxHeadSize',0.8);
end

if ShowBasis
    V   = AMS.facets(index).V;
    adi = V(:,1);
    adj = V(:,2);
    adk = V(:,4);


            
    a = abc(1); b = abc(2); c = abc(3);
    
    % fprintf("a = %d \nb = %d \nc = %d \n",a,b,c);
    %             fprintf("adi = \n")
    %             disp(adi)
    %             fprintf("adj - adi = \n")
    %             disp(adj-adi);
    %             fprintf("adk - adi = \n")
    %             disp(adk - adi);
    
    M = [adi,  b*(adj-adi),  c*(adk - adi)]; 
    M2 = [adi-center' adj adk];
    hold on


    for n = 1:3
        quiver3(center(1),center(2),center(3), M2(1,n), M2(2,n), M2(3,n), ...
            'AutoScale','off','Color',[0,0,n*0.33],'LineWidth',1.5,'MaxHeadSize',0.8);
        % quiver3(adi(1),adi(2),adi(3), M(1,n), M(2,n), M(3,n), ...
        %     'AutoScale','off','Color',[0,0,n*0.33],'LineWidth',1.5,'MaxHeadSize',0.8);
    
        % quiver3(center(1),center(2),center(3), M2(1,n), M2(2,n), M2(3,n), ...
        %     'AutoScale','off','Color',[n*0.3,0,0],'LineWidth',1.5,'MaxHeadSize',0.8);
    end
end





%Grids and axis
box on; grid on
ax = gca;
ax.GridColor = GridColor; ax.GridAlpha = GridAlpha; ax.GridLineStyle = GridStyle;
ax.MinorGridLineStyle = GridStyle; ax.XMinorGrid='on'; ax.YMinorGrid='on'; ax.ZMinorGrid='on';
ax.FontSize = FontSize;

lbl = getfield_with_default(opts,'AxisLabels',{'F_{x} (N)','F_{y} (N)','T_{z} (N·m)'});
xlabel(lbl{1}); ylabel(lbl{2}); zlabel(lbl{3});
title('Attainable Moment Set');

legend([producedHandle,desiredHandle],'moment produced','desired moment')


%Fance lights
if LightingOn
    camlight('headlight'); material dull; lighting gouraud;
end

end

%Options helper
function v = getfield_with_default(s,name,def)
if isfield(s,name), v = s.(name); else, v = def; end
end
