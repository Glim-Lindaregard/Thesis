% cfg = config();
% cfg.A(:,2) = 0;
% U = Copy_of_buildAMS_row(cfg);
% 
% Utest = 0.7*[0 1 1 1 1 1 1 0]';
% index = 0
% 
% for i = 1:size(U,3)
%     Uk = U(:,:,i);                % 8×4 matrix
%     % Check if any column matches Utest exactly
%     if any(all(Uk == Utest, 1))
%         index = i
%         break                           % stop at first match
%     end
% end



%SEE IF BROKEN IS ON
% U = simout.logsout.get('u').Values.Data;
% for i = 1:51
%     u(i) = U(broken,1,i);
% end
% figure;
% plot(1:length(u),u);



% logs = simout.logsout;
% u_ts = logs.get('u');  if isempty(u_ts), u_ts = logs.get('U'); end
% u_ts = u_ts.Values;
% 
% D = u_ts.Data;              % size = [8,1,51]
% if ndims(D)==3 && size(D,2)==1
%     D = squeeze(D);         % -> [8,51]
% end
% if size(D,1)==8             % make time the first dim
%     D = D.';                % -> [51,8]
% end
% 
% tol = 1e-9;
% umin = min(D,[],1);
% umax = max(D,[],1);
% nearZero_always = find(abs(umin)<tol & abs(umax)<tol)
% 
% u_first = D(1,:);
% u_last  = D(end,:);
% disp(table((1:numel(u_first))', u_first(:), u_last(:), ...
%     'VariableNames', {'idx','u_first','u_last'}))
% 


% cfg = config();
% m   = size(cfg.A,2);
% 
% disp('--- THRUSTER GEOMETRY CHECK ---');
% T = table((1:m)', cfg.pos(:,1), cfg.pos(:,2), cfg.beta(:), ...
%     'VariableNames', {'idx','x','y','beta'});
% disp(T);
% 
% figure; hold on; axis equal;
% title('Config geometry (body frame)');
% xlabel('x_b →'); ylabel('y_b ↑');
% for k = 1:m
%     p = cfg.pos(k,:);
%     d = [cos(cfg.beta(k)) sin(cfg.beta(k))];
%     quiver(p(1), p(2), 0.3*d(1), 0.3*d(2), 0, 'LineWidth',1.5, 'Color',[0.3 0.9 0.3]);
%     text(p(1), p(2), sprintf('%d',k), 'VerticalAlignment','bottom','Color',[0 0.9 0]);
% end
% plot(0,0,'ko','MarkerFaceColor','k');
% grid on;


% pull first pose
st = simout.logsout.get('State').Values;
x0 = st.Data(1,1); y0 = st.Data(1,2); th0 = st.Data(1,3);

% config geometry (body frame)
pos  = cfg0.pos;           % 8×2
beta = cfg0.beta(:);       % 8×1

% body->world at the first frame
R = [cos(th0) -sin(th0); sin(th0) cos(th0)];
thr_w = (R*pos.').'+[x0 y0];   % 8×2

% classify quadrants relative to (x0,y0)
q = sign(thr_w - [x0 y0]);      % +1 or -1 per axis
Qpp = find(q(:,1)>0 & q(:,2)>0);   % (+,+)  top-right
Qmp = find(q(:,1)<0 & q(:,2)>0);   % (-,+)  top-left
Qmm = find(q(:,1)<0 & q(:,2)<0);   % (-,-)  bottom-left
Qpm = find(q(:,1)>0 & q(:,2)<0);   % (+,-)  bottom-right

disp(struct('top_right',[Qpp(1),Qpp(2)],'top_left',[Qmp(1), Qmp(2)],'bottom_left',Qmm,'bottom_right',Qpm))
