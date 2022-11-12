function initializeRFPositions(obj)
    % Import mRGC RF positions, for passed eccentricity, size, and
    % eye. Get a little larger retion and crop after we compute the
    % cone spacing below
    [obj.rgcRFpositionsMicrons, obj.rgcRFpositionsDegs] = retinalattice.import.finalMRGCPositions(...
        obj.sourceLatticeSizeDegs, obj.eccentricityDegs, obj.sizeDegs*2.0, obj.whichEye);   
    
    
    if (~isempty(obj.micronsPerDegreeApproximation))
        % Convert positions from degs to microns using the passed
        % microns/deg approximation
        obj.rgcRFpositionsMicrons = obj.rgcRFpositionsDegs * obj.micronsPerDegreeApproximation;
    end

    
    % First crop area that is 10% wider than the desired area
    diff = abs(bsxfun(@minus, obj.rgcRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.55*obj.sizeDegs(1)) & (diff(:,2) <= 0.55*obj.sizeDegs(2)));
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(idx,:);
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(idx,:);
    
    % Compute rf spacings from positions
    obj.rgcRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsMicrons);
    obj.rgcRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsDegs);

    % Crop to desired ROI in degs
    diff = abs(bsxfun(@minus, obj.rgcRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.5*obj.sizeDegs(1)) & (diff(:,2) <= 0.5*obj.sizeDegs(2)));
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(idx,:);
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(idx,:);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(idx);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(idx);
end