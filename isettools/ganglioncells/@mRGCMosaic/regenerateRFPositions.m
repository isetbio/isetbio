function regenerateRFPositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile)

    % Regenerate lattice whose FOV is large enough to encopass the desired size at the desired eccentricity
    fovDegs = sqrt(sum(obj.eccentricityDegs.^2,2)) + max(obj.sizeDegs)*1.3;
    
    obj.rgcRFpositionsMicrons = retinalattice.generatePatch(fovDegs, ...
        'midget ganglion cells', obj.whichEye, exportHistoryToFile, ...
        visualizeConvergence, ...
        obj.useParfor, maxIterations, ...
        'customDegsToMMsConversionFunction', obj.customDegsToMMsConversionFunction, ...
        'customRFspacingFunction', [], ... % do not use a custom mRGC spacing function
        'customMinRFspacing', [] ...       % do not use a custom min mRGC spacing
        );

    % Convert to degs
    obj.rgcRFpositionsDegs = obj.customMMsToDegsConversionFunction(obj.rgcRFpositionsMicrons*1e-3);

    % Compute spacings
    obj.rgcRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsDegs);
    obj.rgcRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsMicrons);
    
    % Crop to desired ROI in degs
    diff = abs(bsxfun(@minus, obj.rgcRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.55*obj.sizeDegs(1)) & (diff(:,2) <= 0.55*obj.sizeDegs(2)));
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(idx,:);
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(idx,:);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(idx);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(idx);
    
end