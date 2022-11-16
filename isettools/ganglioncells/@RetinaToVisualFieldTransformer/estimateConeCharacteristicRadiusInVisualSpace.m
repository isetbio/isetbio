function dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj, theTargetPositionDegs, ...
    simulateCronerKaplanEstimation, neighboringConesNum)

    theConeMosaic = obj.theConeMosaic;
    thePSFData = obj.theVlambdaWeightedPSFData;

    conesNum = numel(theConeMosaic.coneTypes);
    if (conesNum < max([6 neighboringConesNum]))
        fprintf(2, 'The mosaic contain less than the desired number of cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
    
        % Return struct
        dStruct.conesNumInRetinalPatch = 0;
        dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end

    % Find the neighboringConesNum closest (to the target position) cones
    [~,idx] = MosaicConnector.pdist2(theConeMosaic.coneRFpositionsDegs, [], ...
        'fromPosition', theTargetPositionDegs, ...
        'smallest', max([6 neighboringConesNum]) ...
        );
    
    % Estimate mean anatomical cone aperture from the closest (to the target position) cones
    meanConeApertureDegs = mean(theConeMosaic.coneApertureDiametersDegs(idx));
    anatomicalConeCharacteristicRadiusDegs = theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * meanConeApertureDegs;

    hFig = figure(1); clf;
    [visualConeCharacteristicRadiusDegs, bestHorizontalResolutionRotationDegs] = ...
        RetinaToVisualFieldTransformer.analyzeVisuallyProjectedConeAperture(theConeMosaic, ...
            meanConeApertureDegs, thePSFData, simulateCronerKaplanEstimation, hFig);
    
    % Return struct
    dStruct.conesNumInRetinalPatch = conesNum;
    %dStruct.indicesOfConesSortedWithDistanceToTargetRFposition = idx(1:neighboringConesNum);
    dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
    dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
    dStruct.bestHorizontalResolutionRotationDegs = bestHorizontalResolutionRotationDegs;
end
