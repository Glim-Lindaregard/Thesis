cfg = config();
AMS = buildAMS_row(cfg);

Utest = 0.7*[1 1 0 1 0 0 1 0]';
index = 0;

for i = 1:numel(AMS.facets)
    U = AMS.facets(i).U;                % 8×4 matrix
    % Check if any column matches Utest exactly
    if any(all(U == Utest, 1))
        index = i;
        break                           % stop at first match
    end
end
