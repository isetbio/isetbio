function initializeRFPositions(obj)

    % Convert degs to retinal microns
    eccMicrons  = 1e3 * obj.customDegsToMMsConversionFunction(obj.eccentricityDegs); %RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = obj.sizeVisualDegsToSizeRetinalMicrons(obj.sizeDegs*2.0, sqrt(sum(obj.eccentricityDegs.^2,2)));
   
    % Import RF positions from pre-generated lattice
    [obj.rgcRFpositionsMicrons, obj.rgcRFpositionsDegs] = retinalattice.import.finalMRGCPositions(...
        obj.sourceLatticeSizeDegs, eccMicrons, sizeMicrons, obj.whichEye, ...
        obj.customMMsToDegsConversionFunction);  

    % Crop area that is 10% wider than the desired area
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