clear; clc; close all; 


%--- Import model parameters (edit in config_example) ---
cfg = config();


% --- A matrix from Tang paper ---
Atest = [
 -0.18548  -0.18548   0.18548   0.18548   0.31264   0.31264  -0.31264  -0.31264;
  0.31264  -0.31264  -0.31264   0.31264  -0.18548   0.18548   0.18548  -0.18548;
 -0.10173   0.10173  -0.10173   0.10173   0.10173  -0.10173   0.10173  -0.10173
];

%cfg.A = Atest;



ad = [1,1,0]';
[U,A,norms,center] = Copy_of_buildAMS_row(cfg);
[Uout,index,abc] = findUfromAd_DA(ad,U,A);

aProduced = A*Uout;

fprintf("The desired moment was: %d %d %d \n",ad(1),ad(2),ad(3));


fprintf("And the produced moment is: %d %d %d \n",aProduced(1),aProduced(2),aProduced(3));

fprintf("This was done using:\n");
fprintf(" %.2f\n",Uout');


%---Visualize AMS facets---
VisOpts = struct( ...
    'FaceColor', [0.78 1 0.92], ... % pastel cyan
    'EdgeColor', 1*[1 1 1], ...
    'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0.8, ...
    'LineWidth', 0.6, ...
    'BackgroundColor', 0.2*[1 1 1], ...
    'GridColor', 0.4*[0.9 0.9 0.9], ...
    'GridAlpha', 1, ...
    'GridLineStyle', '-', ...
    'UseOctantColors', false, ...
    'Lighting', true, ...
    'ShowNormals', false,...
    'fps', 10, ...
    'Fancy', false, ...
    'Index', index, ...
    'ShowProduced', true,...
    'ShowDesired', true,...
    'ShowBasis', true);


visualizeAMS(U,A,norms,VisOpts,aProduced,ad,abc);

visualizeSlider(cfg,Uout);



