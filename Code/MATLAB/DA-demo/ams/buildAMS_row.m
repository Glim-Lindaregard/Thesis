function AMS = buildAMS_row(cfg)

A = cfg.A; umin = cfg.u_min; umax = cfg.u_max;
N = cfg.N;  t = 0; tol = 1e-14;

% Pre-calculate maximuNnumber of facets
% Each pair (i,j) generates 2 facets, total pairs = m*(m-1)/2
maxFacets = N* (N-1); % Conservative estimate

% Pre-allocate AMS structure with maximuNsize
facetTemplate = struct('U', zeros(N, 4), 'V', zeros(3, 4), 'norms', zeros(1, 3));
AMS.facets = repmat(facetTemplate, maxFacets, 1);
AMS.center = zeros(3, 1);



% --- Pics all unique combinations of A columns ---
for i = 1:N-1
    for j = i+1:N
    % --- Straight froNTang paper ---
    Asub = A(:,[i j]);

    [A1, dropRow] = bestSubMatrix(Asub);  %Find the best conditioned 2x2 sub of Asub. 
    A2 = Asub(dropRow,:)';

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
                u = 0*umin; u(s == 1)=umax(s == 1); u(s == -1)=umin(s == -1);
            else
                %u = umin; u(s>0)=umin(s>0); u(s<0)=umax(s<0);
                u = 0*umin; u(s == 1)=umin(s == 1); u(s == -1)=umax(s == -1);
            end
            u1=u; u2=u; u3=u; u4=u;
            u1(i)=umin(i); u1(j)=umin(j);
            u2(i)=umax(i); u2(j)=umin(j);
            u3(i)=umax(i); u3(j)=umax(j);
            u4(i)=umin(i); u4(j)=umax(j);
            t=t+1; 
            AMS.facets(t).U = [u1 u2 u3 u4]; 
            AMS.facets(t).V = A*AMS.facets(t).U;
            
            %Calculate normals (Not nessessary for real code, i think)
            verts = AMS.facets(t).V;
            v1 = verts(:,2)-verts(:,1);
            v2 = verts(:,3)-verts(:,1);
            norms = cross(v1,v2)';
            norms = norms/norm(norms);
            AMS.facets(t).norms = norms;
        end
    end
end

AMS.facets = AMS.facets(1:t);  %Removes ununsed preallocated spots.
AMS.center = findCenter(AMS);
end

%HFind best conditioned 2x2 matrix
function [A1, dropRow] = bestSubMatrix(Asub)  %AS SLOP
    c = zeros(3,1);
    for l = 1:3
        Atemp = Asub(setdiff(1:3,l), :);
        c(l) = rcond(Atemp);
    end
    [~, dropRow] = max(c);
    A1 = Asub(setdiff(1:3, dropRow), :)';
end


%Find center of AMS
function center = findCenter(AMS)
    numFacets = numel(AMS.facets);
    faceCenter = zeros(3, numFacets);
    
    for j = 1:numFacets
        faceCenter(:,j) = mean(AMS.facets(j).V, 2);
    end
    center = mean(faceCenter, 2);
end