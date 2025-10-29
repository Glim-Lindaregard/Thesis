function [uOut,index,x] = findUfromAd_DA(ad ,U,A)
    found = 0; k = 0;
    tol = 1e-15;
    N = size(A,2);

    % Initialize all outputs at the start
    uOut = zeros(N, 1);
    index = 0;
    x = zeros(3, 1);
    while ~found
        k = k+1;
        Uk = U(:,:,k);
        Vk = A*Uk;

        if k == size(U,3)  + 1
            %fprintf("No such moment possible\n")
            uOut = zeros(N,1);
            found = 1;
            index = 0;
            x = 0 * ones(3,1);
            continue;
        end

        adi = Vk(:,1);
        adj = Vk(:,2);
        adk = Vk(:,4);

  
        M = [ad,  adi - adj,  adi - adk]; 
        
        if rcond(M) < tol
            %fprintf("M singular at %d where the normal is %d %d %d\n",...
            %   k,AMS.facets(k).norms(1),AMS.facets(k).norms(2),AMS.facets(k).norms(3));
            continue;
        else   
            x = M\adi; 
        end
        
        a = x(1); b = x(2); c = x(3);

        if a >= 0 && b >= 0 && c >= 0 && b <= 1 && c <= 1

            ui = U(:,1,k);
            uj = U(:,2,k);
            uk = U(:,4,k);
            uStar = ui + b*(uj-ui) + c*(uk-ui);
            if a >= 1
                uOut = uStar / a;
            else
                uOut = uStar;
            end
            found = 1;
            index = k;
        end
    end
end