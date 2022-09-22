function dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj, theTargetPositionDegs, simulateCronerKaplanEstimation)

    theConeMosaic = obj.theConeMosaic;
    thePSFData = obj.thePSFData;
    coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;

    conesNum = numel(theConeMosaic.coneTypes);
    if (conesNum == 0)
        fprintf(2, 'The mosaic contain no cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
    
        % Return struct
        dStruct.conesNumInRetinalPatch = 0;
        dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end

    % Find the 6 closest (to the target position) cones
    [~,idx] = MosaicConnector.pdist2(theConeMosaic.coneRFpositionsDegs, [], ...
        'fromPosition', theTargetPositionDegs, ...
        'smallest', 6 ...
        );
    
    % Estimate mean anatomical cone aperture from the 6 closest (to the target position) cones
    meanConeApertureDegs = mean(theConeMosaic.coneApertureDiametersDegs(idx));
    anatomicalConeCharacteristicRadiusDegs = coneCharacteristicRadiusConversionFactor * meanConeApertureDegs;

    hFig = figure(1); clf;
    visualConeCharacteristicRadiusDegs = RetinaToVisualFieldTransformer.analyzeVisuallyProjectedConeAperture(...
        anatomicalConeCharacteristicRadiusDegs, thePSFData, simulateCronerKaplanEstimation, hFig);
    
    % Return struct
    dStruct.conesNumInRetinalPatch = conesNum;
    dStruct.indicesOfConesSortedWithDistanceToTargetRFposition = idx;
    dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
    dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
end
