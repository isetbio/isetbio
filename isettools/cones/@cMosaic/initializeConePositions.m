function initializeConePositions(obj)
    % Import cone RF positions, for passed eccentricity, size, and
    % eye. Get a little larger retion and crop after we compute the
    % cone spacing below
    [obj.coneRFpositionsMicrons, obj.coneRFpositionsDegs] = retinalattice.import.finalConePositions(...
        obj.sourceLatticeSizeDegs, obj.eccentricityDegs, obj.sizeDegs*2.0, obj.whichEye);   
    
    
    if (~isempty(obj.micronsPerDegreeApproximation))
        % Convert positions from degs to microns using the passed
        % microns/deg approximation
        obj.coneRFpositionsMicrons = obj.coneRFpositionsDegs * obj.micronsPerDegreeApproximation;
    end

    
    % First crop area that is 10% wider than the desired area
    diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.55*obj.sizeDegs(1)) & (diff(:,2) <= 0.55*obj.sizeDegs(2)));
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);
    
    % Compute cone spacings from positions
    obj.coneRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);
    obj.coneRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsDegs);

    % Crop to desired ROI in degs
    diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.5*obj.sizeDegs(1)) & (diff(:,2) <= 0.5*obj.sizeDegs(2)));
    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
    obj.coneRFspacingsDegs = obj.coneRFspacingsDegs(idx);
    obj.coneRFspacingsMicrons = obj.coneRFspacingsMicrons(idx);
end

