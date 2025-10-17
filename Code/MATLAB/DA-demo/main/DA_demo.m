clear; clc; close all; 


%--- Import model parameters (edit in config_example) ---
cfg = config_example();


% --- A matrix from Tang paper ---
Atest = [
 -0.18548  -0.18548   0.18548   0.18548   0.31264   0.31264  -0.31264  -0.31264;
  0.31264  -0.31264  -0.31264   0.31264  -0.18548   0.18548   0.18548  -0.18548;
 -0.10173   0.10173  -0.10173   0.10173   0.10173  -0.10173   0.10173  -0.10173
];


%--- Construct AMS structure ---
AMSopts = struct('tol_parallell', 1e-12, 'storeIndeces', false);
AMS = buildAMS_null(cfg.A,cfg.u_min,cfg.u_max);
AMS_row = buildAMS(cfg.A, cfg.u_min', cfg.u_max',AMSopts);
%AMS_null = buildAMS_null(Atest ,-cfg.u_max, cfg.u_max);
% --- Normalize AMS verteces ---
%normAMS = normalizeAMS(AMS);

%ad = [1,0,0]';
%ud = findUd(AMS_row,ad)

%fprintf("fx,fy,tau = \n")
%disp(cfg.A * ud);

%---Visualize AMS facets---
VisOpts = struct( ...
    'FaceColor', [0.78 1 0.92], ... % pastel cyan
    'EdgeColor', 0.2*[1 1 1], ...
    'FaceAlpha', 1, ...
    'EdgeAlpha', 0.85, ...
    'LineWidth', 1.2, ...
    'BackgroundColor', 0.2*[1 1 1], ...
    'GridColor', [0.9 0.9 0.9], ...
    'GridAlpha', 1, ...
    'GridLineStyle', '--', ...
    'UseOctantColors', false, ...
    'Lighting', true, ...
    'ShowNormals', false ...
    );
%figure(1)
%visualizeAMS(AMS_row,VisOpts);

%figure(2)
%visualizeSlider(cfg,ud,ad);
%visualizeAMS(AMS_row,VisOpts);



