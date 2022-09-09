function dStruct = estimateConeCharacteristicRadiusInVisualSpace(theConeMosaic, thePSFData, theTargetPositionDegs)
    
    conesNum = numel(theConeMosaic.coneTypes);
    if (conesNum == 0)
        fprintf(2, 'The mosaic contain no cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
            
        % Return struct
        dStruct.conesNumInRetinalPatch = 0;
        dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end

    % Sort cones according to their distance from the mosaic center
    coneDistancesFromTargetPosition = sqrt(sum(bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, theTargetPositionDegs).^2,2));
    [~,idx] = sort(coneDistancesFromTargetPosition, 'ascend');

    % Estimate mean anatomical cone aperture
    conesNumToUse = min([conesNum 6]);
    meanConeApertureDegsInMosaicCenter = mean(theConeMosaic.coneApertureDiametersDegs(idx(1:conesNumToUse)));
    anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * meanConeApertureDegsInMosaicCenter;

    hFig = [];
    videoOBJ = [];
    pdfFileName = '';
    visualConeCharacteristicRadiusDegs = RetinaToVisualFieldTransformer.analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, hFig, videoOBJ, pdfFileName);

   % Return struct
   dStruct.conesNumInRetinalPatch = conesNum;
   dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
   dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
end