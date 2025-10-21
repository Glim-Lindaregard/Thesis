function visualizeAMS(AMS, opts)

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
fps     = getfield_with_default(opts,'fps',10);

figure('Color',BackgroundColor); hold on
count = numel(AMS.facets);

axis manual; axis equal; axis vis3d; ...
    xlim([-2,2]); ylim([-2,2]); zlim([-2,2]);
view(45,30) 

center = zeros(3,1);
%ENABLE IS CONTER OF AMS IS NOT  0 0 0, or set center to known center.
% for j = 1:count
% faceCenter = mean(AMS.facets(j).V');
% center = center + faceCenter;
% end

for k = 1:count
    %Get a facets verteces.
    verts = AMS.facets(k).V;
    
    %extract verteces
    A = verts(:,1)'; B = verts(:,2)'; C = verts(:,3)'; D = verts(:,4)';
    
    %tris = [A; B; C; A; C; D];  If plotting looks weird, might be this.

    tris = [A; B; C; D];
    
    %f = [1 2 3; 4 5 6];   And this

    f = [1 2 3 4];
    
    patch('Faces',f,'Vertices',tris, ...
          'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
          'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',LineWidth, ...
          'EdgeLighting','none','FaceLighting','flat');
    
    %Fancy camera spin
    for i = 1:fps
    camorbit(360/(fps*count),0,'data'); 
    drawnow
    %pause(0.002)
    end
    drawnow limitrate

    %Edges %Put back on if edges look wierd
    % outline = [A;B;C;D;A];
    % plot3(outline(:,1),outline(:,2),outline(:,3), ...
    %       'Color',EdgeColor,'LineWidth',LineWidth,'LineJoin','round');

    %Show normals
    if ShowNormals
        n = AMS.facets(k).norms;
        ctr = mean(verts,2);
        %Make sure normals point in "correct" direction
        if dot(ctr-center,n)<0
            n = -n;
        end
        quiver3(ctr(1),ctr(2),ctr(3),0.2*n(1),0.2*n(2),0.2*n(3), ...
                    'AutoScale','off','Color','r','LineWidth',1.5,'MaxHeadSize',0.8);
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

%Fance lights
if LightingOn
    camlight('headlight'); material dull; lighting gouraud;
end

end

%Options helper
function v = getfield_with_default(s,name,def)
if isfield(s,name), v = s.(name); else, v = def; end
end
