function limitTestAMS(cfg)
    failed = 0;
    ads = [                 % 1. pure zero vector
    1     0     0;      % 2. +Fx only
   -1     0     0;      % 3. -Fx only
    0     1     0;      % 4. +Fy only
    0    -1     0;      % 5. -Fy only
    0     0     1;      % 6. +tau only
    0     0    -1;      % 7. -tau only
    1     1     0;      % 8. diagonal force (45°)
   -1    -1     0;      % 9. opposite diagonal
    1     0     1;      % 10. +Fx and +tau
    0     1    -1;      % 11. +Fy and -tau
    1    -1     0;      % 12. opposite force components
    0.7   0.7   0.7;    % 13. all positive mix
   -0.7  -0.7   0.7;    % 14. mixed signs
    2     0     0;      % 15. beyond nominal Fx range
    0     2     0;      % 16. beyond nominal Fy range
    0     0     2;      % 17. beyond nominal torque range
    0.01  0.01  0;      % 18. near-zero but nonzero
    1e-8  1e-8  1e-8;   % 19. effectively zero (numerical tolerance test)
    ];
    
    
    for i = 1:length(ads(:,1))
        %--- Construct AMS structure ---
        ad = ads(i,:)';
        
        m = cfg.m;
        AMS = buildAMS_row(cfg);
    
        if numel(AMS.facets) ~= m^2 -m
            fprintf("LIMIT TESTER: Some facets failed to ... " + ...
                "build, %d facets where built\", nume(AMS.facets));
        end
    
    
        [U,index] = findUfromAd_DA(ad, AMS);
        
        if ~index
            fprintf("LIMIT TESTER: Test nr:%d has failed\n",index);
            failed = failed + 1;
        end
        
        % aProduced = cfg.A*U;
        % 
        % fprintf("The desired moment was: %d %d %d \n",ad(1),ad(2),ad(3));
        % 
        % 
        % fprintf("And the produced moment is: %d %d %d \n",aProduced(1),aProduced(2),aProduced(3));
        % 
        % fprintf("This was done using:\n");
        % fprintf(" %.1f\n",U');
    end
    fprintf("%d moment(s) were not attained: TEST PASSED\n",failed);
end