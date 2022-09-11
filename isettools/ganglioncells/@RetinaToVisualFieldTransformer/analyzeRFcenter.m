% Method to compute the cone map for the RF center and its corresponding Gaussian characteristic radius
function [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs, ...
    visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
    visualRFcenterXYpos, visualRFcenterOrientationDegs] = analyzeRFcenter(obj, ...
        indicesOfConesPooledByTheRFcenter, ...
        weightsOfConesPooledByTheRFcenter, ...
        spatialSupportDegs)
    
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

    theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(spatialSupportDegs(:,1), spatialSupportDegs(:,2), visualRFcenterConeMap, ...
        'flatTopGaussian', obj.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation, ...
        'forcedOrientationDegs', forcedOrientationDegs, ...
        'globalSearch', false);

    % Return the characteristic radius in degrees
    % The sum of the 2 radii
    %visualRFcenterCharacteristicRadiusDegs = sqrt(sum(theFittedGaussian.characteristicRadii.^2,2))/sqrt(2);

    % The min of the 2 radii
    %visualRFcenterCharacteristicRadiusDegs = min(theFittedGaussian.characteristicRadii)
    %visualRFcenterCharacteristicRadiusDegs = mean(theFittedGaussian.characteristicRadii)
    visualRFcenterCharacteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));
    
    visualRFcenterCharacteristicRadiiDegs = theFittedGaussian.characteristicRadii;
    visualRFcenterFlatTopExponents = theFittedGaussian.flatTopExponents;
    visualRFcenterXYpos = theFittedGaussian.xyCenter;
    visualRFcenterOrientationDegs = theFittedGaussian.rotationDegs;
   
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

