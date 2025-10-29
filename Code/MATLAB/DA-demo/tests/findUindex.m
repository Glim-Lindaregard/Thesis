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
U = simout.logsout.get('u').Values.Data;
for i = 1:51
    u(i) = U(broken,1,i);
end
figure;
plot(1:length(u),u);
