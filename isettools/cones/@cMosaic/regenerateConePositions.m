function regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile)
    
    if ((obj.eccentricityDegs(1) ~= 0) || (obj.eccentricityDegs(2) ~= 0))
        error('When computing mesh from scratch, ''eccentricityDegs'' must be set to [0 0].');
    end
    
    % Regenerate lattice over a somewhat larger region, centered at (0,0)
    obj.coneRFpositionsMicrons = retinalattice.generatePatch(max(obj.sizeDegs)*1.3, ...
        'cones', obj.whichEye, exportHistoryToFile, visualizeConvergence, obj.useParfor, maxIterations);
    
    % Convert to degs
    obj.coneRFpositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(obj.coneRFpositionsMicrons*1e-3);
    
    % Compute spacings (which determine apertures)
    obj.coneRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsDegs);
    obj.coneRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);
    
    % Crop to desired ROI in degs
    diff = abs(obj.coneRFpositionsDegs);
    idx = find((diff(:,1) <= 0.5*obj.sizeDegs(1)) & (diff(:,2) <= 0.5*obj.sizeDegs(2)));
    
    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
    obj.coneRFspacingsDegs = obj.coneRFspacingsDegs(idx);
    obj.coneRFspacingsMicrons = obj.coneRFspacingsMicrons(idx);
    
end
