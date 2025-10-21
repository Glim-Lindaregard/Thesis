function AMS = buildAMS_row(cfg)
%Assign cfg variables etc
A = cfg.A; umin = cfg.u_min; umax = cfg.u_max;
[~,m] = size(A); t = 0;
tol = 1e-12;

% --- Pics all unique combinations of A columns ---
for i = 1:m-1
    for j = i+1:m
    % --- Straight from Tang paper ---
    Asub = A(:,[i j]); 
    [A1, dropRow] = bestSubMatrix(Asub);  %Find the best conditioned 2x2 sub of Asub. 
    A2 = Asub(dropRow,:)';

    n = zeros(3,1);
    n(setdiff(1:3,dropRow)) = -A1 \ A2;
    n(dropRow) = 1;

    if norm(n) < tol
        % columns nearly parallel; skip this pair
        frpintf("Cross product too small");
        continue
    end

    % --- For each pair, generate facet vertices ---
    s = sign(n.'*A); 
    s([i j]) = 0;

        for which = 1:2
            %Calculate verteces
            if which == 1
                u = umin; u(s>0)=umax(s>0); u(s<0)=umin(s<0);
            else
                u = umin; u(s>0)=umin(s>0); u(s<0)=umax(s<0);
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
end

%Helper to find best conditioned 2x2 matrix
function [A1, dropRow] = bestSubMatrix(Asub)
    c = zeros(3,1);
    for l = 1:3
        Atemp = Asub(setdiff(1:3,l), :);
        c(l) = rcond(Atemp);
    end
    [~, dropRow] = max(c);
    A1 = Asub(setdiff(1:3, dropRow), :)';
end