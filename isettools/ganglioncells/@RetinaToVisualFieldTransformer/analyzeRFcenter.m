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
    [retinalRFcenterConeMap, anatomicalRFcenterCharacteristicRadiusDegs] = computeRetinalRFcenterMapFromInputConesInArcMinSupport(obj, ...
        indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs, ...
        flatTopGaussian, forcedOrientationDegs);
    
    
    % Compute the visual RF center cone map via convolution of the retinalRFcenterConeMap with the PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.theVlambdaWeightedPSFData.vLambdaWeightedData, 'same');

   
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

        % Fit a 1D Gaussian line weighting function to the 1D profile 
        % (integration along the Y-dimension of the 2D visually projected
        % cone aperture map)
        visualRFcenterConeMapProfile = sum(rotatedVisualRFcenterConeMap,1);
        theFittedGaussianLineWeightingFunction = RetinaToVisualFieldTransformer.fitGaussianLineWeightingFunction(...
            obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs, visualRFcenterConeMapProfile);

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
            'globalSearch', false);

        % The sqrt(product) of the 2 radii
        visualRFcenterCharacteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));
    
        visualRFcenterCharacteristicRadiiDegs = theFittedGaussian.characteristicRadii;
        visualRFcenterFlatTopExponents = theFittedGaussian.flatTopExponents;
        visualRFcenterXYpos = theFittedGaussian.xyCenter;
        visualRFcenterOrientationDegs = theFittedGaussian.rotationDegs;
    end

   
end

function [retinalRFcenterConeMap, anatomicalRFcenterCharacteristicRadiusDegs] = computeRetinalRFcenterMapFromInputConesInArcMinSupport(obj, ...
    theConeIndices, theConeWeights, spatialSupportDegs, flatTopGaussian, forcedOrientationDegs)

    % Compute the retinal RF center cone map
    theConeCharacteristicRadiiDegs = obj.coneCharacteristicRadiusConversionFactor * obj.theConeMosaic.coneApertureDiametersDegs(theConeIndices);
    theConePositionsDegs = obj.theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);
    retinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConePositionsDegs, theConeWeights, spatialSupportDegs);

    theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
        spatialSupportDegs(:,1), spatialSupportDegs(:,2), retinalRFcenterConeMap, ...
        'flatTopGaussian', flatTopGaussian, ...
        'forcedOrientationDegs', forcedOrientationDegs, ...
        'globalSearch', false);

    % The sqrt(product) of the 2 radii
    anatomicalRFcenterCharacteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));
end



