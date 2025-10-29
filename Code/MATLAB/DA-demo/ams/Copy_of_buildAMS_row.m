function [U,A,norms,center] = Copy_of_buildAMS_row(cfg)

A = cfg.A; umin = cfg.u_min; umax = cfg.u_max;
N = cfg.N;  t = 0; tol = 1e-14;

% Pre-calculate maximuNnumber of facets
% Each pair (i,j) generates 2 facets, total pairs = m*(m-1)/2
maxFacets = N* (N-1); % Conservative estimate

% Pre-allocate AMS structure with maximuNsize
U = zeros(N,4,maxFacets);
V = zeros(3,4,maxFacets);
norms = zeros(1,3,maxFacets);
%center = zeros(3, 1);


zeroIndex = find(all(abs(A) < 0.01),1);

if isempty(zeroIndex)
    zeroIndex = -1;
end

% --- Pics all unique combinations of A columns ---
for i = 1:N-1
    if i == zeroIndex
        continue
    end
    for j = i+1:N
        if j == zeroIndex
            continue
        end
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
            U(:,:,t) = [u1 u2 u3 u4]; 
            V(:,:,t) = A*U(:,:,t);
            
            %Calculate normals (Not nessessary for real code, i think)
            verts = V(:,:,t);
            v1 = verts(:,2)-verts(:,1);
            v2 = verts(:,3)-verts(:,1);
            normal = cross(v1,v2)';
            normal = normal/norm(normal);
            norms(:,:,t) = normal;
        end
    end
end

V = V(:,:,1:t);  %Removes ununsed preallocated spots.
center = findCenter(V);
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
function center = findCenter(V)
    numFacets = size(V,3);
    faceCenter = zeros(3, numFacets);
    
    for j = 1:numFacets
        faceCenter(:,j) = mean(V(:,:,j), 2);
    end
    center = mean(faceCenter, 2);
end