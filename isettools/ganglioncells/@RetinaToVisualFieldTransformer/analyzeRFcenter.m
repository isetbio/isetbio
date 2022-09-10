function [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs] = analyzeRFcenter(obj, ...
    indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs)
    
    % Compute the visual RF center map
    visualRFcenterConeMap = computeVisualRFcenterMapFromInputConesInArcMinSupport(obj, ...
        indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs);
    
    % Compute its characteristic radius by fitting a flat top gaussian
    % ellipsoid to it.
    if (numel(indicesOfConesPooledByTheRFcenter) == 2)
       % We are using a circularly summetric PSF, so force the
       % orientation to match the orientation of the 2 cones
       cone1RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter(1),:);
       cone2RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter(2),:);
       deltaY = cone2RFpos(2)-cone1RFpos(2);
       deltaX = cone2RFpos(1)-cone1RFpos(1);
       forcedOrientationDegs = -atan2d(deltaY, deltaX);
    else
       forcedOrientationDegs = [];
    end

    flatTopGaussian = true;
    theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(spatialSupportDegs(:,1), spatialSupportDegs(:,2), visualRFcenterConeMap, ...
        'flatTopGaussian', flatTopGaussian, ...
        'forcedOrientationDegs', forcedOrientationDegs, ...
        'globalSearch', false);

    % Return the characteristic radius in degrees
    visualRFcenterCharacteristicRadiusDegs = sqrt(sum(theFittedGaussian.characteristicRadii.^2,2))/sqrt(2);
end

function visualRFcenterConeMap = computeVisualRFcenterMapFromInputConesInArcMinSupport(obj, ...
    theConeIndices, theConeWeights, spatialSupportDegs)

    % Compute the retinal RF center cone map
    theConeCharacteristicRadiiDegs = obj.coneCharacteristicRadiusConversionFactor * obj.theConeMosaic.coneApertureDiametersDegs(theConeIndices);
    theConePositionsDegs = obj.theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);
    retinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConePositionsDegs, theConeWeights, spatialSupportDegs);

    % Compute the visual RF center cone map via convolution of the retinalRFcenterConeMap with the PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.theCircularPSFData.data, 'same');

    % Unit amplitude
    visualRFcenterConeMap = visualRFcenterConeMap ./ max(visualRFcenterConeMap(:));
end

