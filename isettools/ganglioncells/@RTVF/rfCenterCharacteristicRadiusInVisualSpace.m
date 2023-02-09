function visualRFcenterRcDegs = rfCenterCharacteristicRadiusInVisualSpace(obj)

    theConeIndices = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;
    theConeWeights = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;
    theConeApertureAreas = obj.coneMosaic.computeApertureAreasMetersSquared(theConeIndices);
    theEffectiveOSlengthAttenuationFactors = obj.coneMosaic.computeEffectiveOSlengthAttenuationFactors(theConeIndices);
    
    % Re-center the input cone positions
    theConePositionsDegs = obj.coneMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);

    % Compute the retinal cone input map for the RF center
    retinalRFcenterConeMap = obj.retinalSubregionConeMapFromPooledConeInputs(...
        theConePositionsDegs, ...
        theConeApertureAreas, ...
        theEffectiveOSlengthAttenuationFactors, ...
        theConeWeights);

    % Convolve with the LM-weighted PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.spectrallyWeightedPSFData.LMconeWeighted, 'same');

    % Since RF parameters by Croner&Kaplan were based on optimal
    % orientation gratings (highest resolving SF), we first rotate the
    % visualRFcenterConeMap so that it is narrowest along the x-axis,
    % then sum along the y-dimension, and finally fit a line-weighting function 
    % for a Gaussian 

    % First rotate
    debugRadonTransform = false;
    [visualRFcenterConeMapRotated, rotationDegs] = ...
        RTVF.bestHorizontalResolutionRFmap(visualRFcenterConeMap, [], debugRadonTransform);

    % Then integrate along y of the 2-D RF
    visualRFcenterProfileX = sum(visualRFcenterConeMapRotated,1);

    % Finally, fit a 1D Gaussian line weighting function to the 1D profile 
    spatialSupportXdegs = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
    theFittedGaussianLineWeightingFunction = ...
        RTVF.fitGaussianLineWeightingFunction(spatialSupportXdegs, visualRFcenterProfileX);
    
    % Return the visualRF center Rc (degs)
    visualRFcenterRcDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;


    debugFitting = true;
    if (debugFitting)
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 950 950])
        ax = subplot(2,2,1);
        imagesc(...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
            retinalRFcenterConeMap);
        axis 'xy'; axis 'image'
        colormap(gray(1024))
        title('retinal rf center');
    
        ax = subplot(2,2,2);
        imagesc(...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
            visualRFcenterConeMap);
        axis 'xy'; axis 'image'
        title('visual rf center');

        ax = subplot(2,2,3);
        imagesc(...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
            obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
            visualRFcenterConeMapRotated);
        axis 'xy'; axis 'image'
        title(sprintf('visual rf center (rot = %2.1f degs)', rotationDegs));
    
        ax = subplot(2,2,4);
        plot(obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), visualRFcenterProfileX, 'ks');
        hold on;
        plot(spatialSupportXdegs+ meanConePositionDegs(1), theFittedGaussianLineWeightingFunction.profile, 'r-');
        plot([0 visualRFcenterRcDegs]+ meanConePositionDegs(1), max(theFittedGaussianLineWeightingFunction.profile)*exp(-1)*[1 1], 'b--');
        axis 'square'
    end
    

end 

