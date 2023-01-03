% Method to compute the cone map for the RF center and its corresponding Gaussian characteristic radius
function [visualRFcenterCharacteristicRadiusDegs, visualRFcenterConeMap, ...
    visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
    visualRFcenterXYpos, visualRFcenterOrientationDegs, ...
    anatomicalRFcenterCharacteristicRadiusDegs] = analyzeRFcenter(obj, ...
        indicesOfConesPooledByTheRFcenter, ...
        weightsOfConesPooledByTheRFcenter, ...
        spatialSupportDegs)
    
    if (numel(indicesOfConesPooledByTheRFcenter) == 2)
           cone1RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter(1),:);
           cone2RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter(2),:);
           deltaY = cone2RFpos(2)-cone1RFpos(2);
           deltaX = cone2RFpos(1)-cone1RFpos(1);
           forcedOrientationDegs = -atan2d(deltaY, deltaX);
        else
           forcedOrientationDegs = [];
    end

    if (numel(indicesOfConesPooledByTheRFcenter) == 1)
        flatTopGaussian = false;
    else
        flatTopGaussian = obj.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation;
    end

    % Compute the visual RF center map
    [retinalRFcenterConeMap, anatomicalRFcenterCharacteristicRadiusDegs] = computeRetinalRFcenterMapFromInputCones(obj, ...
        indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs, ...
        flatTopGaussian, forcedOrientationDegs);
    
    
    % Compute the visual RF center cone map via convolution of the retinalRFcenterConeMap with the PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.theVlambdaWeightedPSFData.vLambdaWeightedData, 'same');

    % Normalize to unit amplitude
    %visualRFcenterConeMap = visualRFcenterConeMap / max(visualRFcenterConeMap(:));

    if (obj.simulateCronerKaplanEstimation)
        % Since RF parameters by Croner&Kaplan were based on gratings, to
        % approximate this (and to include features of this estimation) we sum
        % along the y-dimension of the visualluy projected cone aperture map
        % and subsequently fit a line-weighting function for a Gaussian 

        % Rotate the targetVisualRFmap so as to maximize horizontal resolution and 
        % retrieve the rotation that maximizes horizontal resolution
        bestHorizontalResolutionRotationDegs = [];
        [rotatedVisualRFcenterConeMap, bestHorizontalResolutionRotationDegs] = ...
            RetinaToVisualFieldTransformer.bestHorizontalResolutionRFmap(visualRFcenterConeMap, bestHorizontalResolutionRotationDegs);

        % Integrate over Y
        visualRFcenterConeMapProfile = sum(rotatedVisualRFcenterConeMap,1);

        % Normalize to unit amplitude
        visualRFcenterConeMapProfile = visualRFcenterConeMapProfile / max(visualRFcenterConeMapProfile(:));

        % Fit a 1D Gaussian line weighting function to the 1D profile 
        % (integration along the Y-dimension of the 2D visually projected
        % cone aperture map)
        theFittedGaussianLineWeightingFunction = RetinaToVisualFieldTransformer.fitGaussianLineWeightingFunction(...
            obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, visualRFcenterConeMapProfile);

        % Compute the corresponding STFs
        [spatialFrequencySupport, visualRFcenterSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, visualRFcenterConeMapProfile);
        [~, visualRFcenterModelSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, theFittedGaussianLineWeightingFunction.profile);
        [~, visualRFcenterModelSTFNormalized] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, theFittedGaussianLineWeightingFunction.profile/max(theFittedGaussianLineWeightingFunction.profile));



        hFig = figure(444); clf;
        set(hFig, 'Position', [10 10 950 1150], 'Name', 'RF center analysis');
        subplot(2,2,1)
        imagesc(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, ...
                obj.theVlambdaWeightedPSFData.spatialSupportForRFmapYdegs, ...
                visualRFcenterConeMap);
        axis 'image';
        set(gca, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'FontSize', 16);
        colormap(gray(1024))
        title('RF center cone map');
        
        subplot(2,2,2);
        imagesc(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, ...
                obj.theVlambdaWeightedPSFData.spatialSupportForRFmapYdegs, ...
                rotatedVisualRFcenterConeMap);
        axis 'image'
        set(gca, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'FontSize', 16);
        colormap(gray(1024))
        title('RF center cone map (rotated)');


        subplot(2,2,3);
        plot(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, visualRFcenterConeMapProfile, 'ko', 'MarkerSize', 12);
        hold on
        plot(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, theFittedGaussianLineWeightingFunction.profile, 'r-', 'LineWidth', 1.5);
        plot(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, theFittedGaussianLineWeightingFunction.profile/max(theFittedGaussianLineWeightingFunction.profile), 'r--', 'LineWidth', 1.0);
        
        legend('RFcenter profile', 'Gaussian fit', 'Gaussian fit (normalized)', 'Location', 'NorthOutside')
        set(gca, 'XLim', 0.1*[-1 1], 'FontSize', 16);
        axis 'square'

        % Plot the STFs
        subplot(2,2,4);
        plot(spatialFrequencySupport, visualRFcenterSTF, 'ko', 'MarkerSize', 12); hold on;
        plot(spatialFrequencySupport, visualRFcenterModelSTF, 'r-', 'LineWidth', 1.5);
        plot(spatialFrequencySupport, visualRFcenterModelSTFNormalized, 'r--', 'LineWidth', 1.0);
        legend('RFcenter profile', 'Gaussian fit', 'Gaussian fit (normalized)', 'Location', 'NorthOutside')
        set(gca, 'XLim', [0.1 100], 'XScale', 'log', 'FontSize', 16);
        axis 'square'

        drawnow;
        
        % Return the characteristic radius in degrees
        visualRFcenterCharacteristicRadiusDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;

        visualRFcenterCharacteristicRadiiDegs = [];
        visualRFcenterFlatTopExponents = [];
        visualRFcenterXYpos = [theFittedGaussianLineWeightingFunction.xo 0];
        visualRFcenterOrientationDegs = bestHorizontalResolutionRotationDegs;

    else

        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportDegs(:,1), spatialSupportDegs(:,2), ...
            visualRFcenterConeMap, ...
            'flatTopGaussian', flatTopGaussian, ...
            'forcedOrientationDegs', forcedOrientationDegs, ...
            'globalSearch', true, ...
            'multiStartsNum', 16);

        % The sqrt(product) of the 2 radii
        visualRFcenterCharacteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));
    
        visualRFcenterCharacteristicRadiiDegs = theFittedGaussian.characteristicRadii;
        visualRFcenterFlatTopExponents = theFittedGaussian.flatTopExponents;
        visualRFcenterXYpos = theFittedGaussian.xyCenter;
        visualRFcenterOrientationDegs = theFittedGaussian.rotationDegs;
    end

end

function [retinalRFcenterConeMap, anatomicalRFcenterCharacteristicRadiusDegs] = ...
    computeRetinalRFcenterMapFromInputCones(obj, ...
    theConeIndices, theConeWeights, spatialSupportDegs, flatTopGaussian, forcedOrientationDegs)

    % Compute the retinal RF center cone map
    theConePositionsDegs = obj.theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);
    retinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        obj.theConeMosaic, ...
        theConePositionsDegs, ...
        theConeWeights, ...
        spatialSupportDegs);

    theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
        spatialSupportDegs(:,1), spatialSupportDegs(:,2), retinalRFcenterConeMap, ...
        'flatTopGaussian', flatTopGaussian, ...
        'forcedOrientationDegs', forcedOrientationDegs, ...
        'globalSearch', false);

    % The sqrt(product) of the 2 radii
    % anatomicalRFcenterCharacteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));

    % The minimum of the 2 radii
    anatomicalRFcenterCharacteristicRadiusDegs = min(theFittedGaussian.characteristicRadii);
end



