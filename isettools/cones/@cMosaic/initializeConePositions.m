function initializeConePositions(obj)
    % Import cone RF positions, for passed eccentricity, size, and
    % eye. Get a little larger retion and crop after we compute the
    % cone spacing below
    [obj.coneRFpositionsMicrons, obj.coneRFpositionsDegs] = retinalattice.import.finalConePositions(...
        obj.sourceLatticeSizeDegs, obj.eccentricityDegs*1.1, obj.sizeDegs*1.3, obj.whichEye);   

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

