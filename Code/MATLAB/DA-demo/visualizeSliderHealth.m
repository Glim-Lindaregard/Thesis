function visualizeSliderHealth(health, ax)
% visualizeSliderHealth(health, ax)
% health: 1x8 or 8x1 with values:
%   1 = healthy   (gray)
%   0 = passive   (red)
%   2 = active    (blue)
%
% Layout per your sketch. Arrows point from outside into the square.

health = health(:).';
if numel(health) ~= 8
    error('health must have 8 elements.');
end

% --- geometry (drawing scale) ---
a = 220;
corners = [ a  a;   -a  a;   -a -a;    a -a];  % [TR; TL; BL; BR]

thr_corner = [1 1 2 2 3 3 4 4];
dir = [ ...
   -1  0;   % T1
    0 -1;   % T2
    0 -1;   % T3
    1  0;   % T4
    1  0;   % T5
    0  1;   % T6
    0  1;   % T7
   -1  0];  % T8

dir = -dir;


pos = corners(thr_corner,:);

% --- plotting frame ---
if nargin < 2 || isempty(ax)
    fig = figure('Color','w');
    ax = axes('Parent', fig);
end

cla(ax);
hold(ax,'on');
axis(ax,'equal');

fig = ancestor(ax,'figure');
set(fig,'Color','w');
set(ax,'Color','w');
                 % axes background
set(ax,'XColor','k','YColor','k');    % axis numbers black

marg = 0.6*a;
xlim(ax, [-a-marg, a+marg]);
ylim(ax, [-a-marg, a+marg]);

rectangle(ax,'Position',[-a -a 2*a 2*a],'LineWidth',2);

% Colors
colHealthy = [0.35 0.35 0.35];
colPassive = [1.00 0.20 0.20];
colActive  = [0.20 0.45 1.00];

% Style knobs
L        = 0.35*a;    % arrow length (from outside toward corner)
ms       = 9;         % marker size
headGap  = 0.12*a;    % pull arrowhead back from corner (MUST be positive)
arrowLW  = 2.5;       % arrow line width
arrowHead= 3.2;       % arrow head size

% Label offsets (bigger + tuned)
labelOffset = 0.30*a * [ ...
     1.2,  0.2;   % T1  (TR, push right)
     0.2,  1.2;   % T2  (TR, push up)
    -0.2,  1.2;   % T3  (TL, push up)
    -1.2,  0.2;   % T4  (TL, push left)
    -1.2, -0.2;   % T5  (BL, push left)
    -0.2, -1.2;   % T6  (BL, push down)
     0.2, -1.2;   % T7  (BR, push down)
     1.2, -0.2];  % T8  (BR, push right)

for i = 1:8
    if health(i) == 0
        c = colPassive;
    elseif health(i) == 2
        c = colActive;
    else
        c = colHealthy;
    end

    % --- Arrow with HEAD on the OUTSIDE ---
    % --- Arrow (HEAD OUTSIDE), full length, not hidden by marker ---
    d = dir(i,:) / max(norm(dir(i,:)), eps);
    
    % Start just outside the marker, at the corner
    tail = pos(i,:) + (0.55*ms) * d;   % push tail away from corner by ~marker radius
    
    % Full arrow length outward
    dx = L * d(1);
    dy = L * d(2);
    
    quiver(ax, tail(1), tail(2), dx, dy, 0, ...
        'LineWidth', arrowLW, 'MaxHeadSize', arrowHead, 'Color', c);



    % --- Marker (draw on top) ---
    plot(ax, pos(i,1), pos(i,2), 'ko', 'MarkerFaceColor', c, 'MarkerSize', ms);

    % --- Label (black, with white background; aligned away from center) ---
    ptxt = pos(i,:) + labelOffset(i,:);

    if ptxt(1) > 0, ha = 'left';  else, ha = 'right'; end
    if ptxt(2) > 0, va = 'bottom'; else, va = 'top'; end

    text(ax, ptxt(1), ptxt(2), sprintf('%d', i), ...
        'FontSize', 14, 'FontWeight','bold', 'Color','k', ...
        'HorizontalAlignment', ha, 'VerticalAlignment', va, ...
        'BackgroundColor','w', 'Margin', 1);
end

% --- Legend: use real (finite) arrows so glyphs show up ---
Lleg = 0.15*a;
hH = quiver(ax, -10*a, -10*a, Lleg, 0, 0, 'Color', colHealthy, 'LineWidth', 2, 'MaxHeadSize', 1.5);
hA = quiver(ax, -10*a, -10*a, Lleg, 0, 0, 'Color', colActive,  'LineWidth', 2, 'MaxHeadSize', 1.5);
hP = quiver(ax, -10*a, -10*a, Lleg, 0, 0, 'Color', colPassive, 'LineWidth', 2, 'MaxHeadSize', 1.5);

lgd = legend(ax, [hH hA hP], {'Healthy','Active','Passive'}, 'Location', 'best');

set(lgd,'TextColor','k','Color','w','EdgeColor',[0.2 0.2 0.2]);

title(ax, 'Thruster health (arrows point into the square)', 'Color','k');
xlabel(ax, 'x', 'Color','k');
ylabel(ax, 'y', 'Color','k');
grid(ax, 'on');

hold(ax,'off');
end
