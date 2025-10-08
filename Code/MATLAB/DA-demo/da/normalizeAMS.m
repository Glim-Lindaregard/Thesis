function normalizedAMS = normalizeAMS(AMS)
    
    numberOfFacets = AMS.num_facets_built;
    normalizedAMS = AMS;

    for i = 1:numberOfFacets

        verts = AMS.facets(i).Uverts;
        tau = AMS.facets(i).Tau;

        % --- Normalize each column vector ---
        nVert = vecnorm(verts);
        nTau = vecnorm(tau);
        nVert(nVert==0) = 1;   %If n = 0 set n = 1 and do nothing with the column.
        nTau(nTau==0) = 1;   %If n = 0 set n = 1 and do nothing with the column.
        normalizedAMS.facets(i).Uverts = verts ./ nVert;
        normalizedAMS.facets(i).Tau = tau ./ nTau;
    end
    
end