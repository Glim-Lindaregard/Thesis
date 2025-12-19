function U = buildAMS_row(cfg)

A = cfg.A; umin = cfg.u_min; umax = cfg.u_max;
N = cfg.N_thrusters;  t = 0; tol = 1e-15;


maxFacets = N* (N-1); % Conservative estimate

% Pre-allocate AMS structure with maximuNsize
U = zeros(N,4,maxFacets);

zeroIndex = find(all(abs(A) < tol));


%If only 1 or 2 thrusters are healthy set AMS to zero. (give up)
if numel(zeroIndex) >= size(A,2) - 1
    U = zeros(8,4,56);
    return;
end

if isempty(zeroIndex)
    zeroIndex = -1;
end

% --- Pics all unique combinations of A columns ---
for i = 1:N-1
    if any(zeroIndex == i)
        continue
    end
    for j = i+1:N
        if any(zeroIndex == j)
            continue
        end

        %Pick Asub = [Ai,Aj] (3x2)
        Asub = A(:,[i j]);
    
    
        [A1, dropRow] = bestSubMatrix(Asub);  %Find the best conditioned 2x2 sub of Asub. 
        A2 = Asub(dropRow,:)';
    
        %Create normal vector
        n = zeros(3,1);
        n(setdiff(1:3,dropRow)) = -A1 \ A2;
        n(dropRow) = 1;

        
        n(abs(n) < tol) = 0;
        if norm(n) < tol
            % columns nearly parallel; skip this pair
            fprintf("Cross product too small");
            continue
        end
    
        % --- For each pair, generate facet vertices ---
        ss = (n.'*A);
        s = zeros(1,N);
        for l = 1:N
            if (ss(l) < tol) && (ss(l) > -tol)
                s(l) = 0;
            elseif ss(l) < -tol
                s(l) = -1;
            else 
                s(l) = 1;
            end
        end
        s([i j]) = 0;

        for which = 1:2
            %Calculate verteces
            if which == 1
                %u = umin; u(s>0)=umax(s>0); u(s<0)=umin(s<0);
                u = umin; u(s == 1)=umax(s == 1); u(s == -1)=umin(s == -1);
            else
                %u = umin; u(s>0)=umin(s>0); u(s<0)=umax(s<0);
                u = umin; u(s == 1)=umin(s == 1); u(s == -1)=umax(s == -1);
            end
            u1=u; u2=u; u3=u; u4=u;
            u1(i)=umin(i); u1(j)=umin(j);
            u2(i)=umax(i); u2(j)=umin(j);
            u3(i)=umax(i); u3(j)=umax(j);
            u4(i)=umin(i); u4(j)=umax(j);
            t=t+1; 
            U(:,:,t) = [u1 u2 u3 u4]; 
        end
    end
end

end

%Find best conditioned 2x2 matrix
function [A1, dropRow] = bestSubMatrix(Asub)
    c = zeros(3,1);
    for l = 1:3
        Atemp = Asub(setdiff(1:3,l), :);
        c(l) = rcond(Atemp);
    end
    [~, dropRow] = max(c);
    A1 = Asub(setdiff(1:3, dropRow), :)';
end