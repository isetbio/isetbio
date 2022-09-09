function [retinalRFparams, weightsComputeFunctionHandle, ...
          targetVisualRFmap, spatialSupportDegs, modelConstants ] = retinalRFparamsForTargetVisualRF(...
                obj,indicesOfConesPooledByTheTargetRFcenter, targetVisualRFDoGparams)
    

    spatialSupportDegs = [obj.thePSFData.supportXdegs(:) obj.thePSFData.supportYdegs(:)];

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs] = ...
       analyzeRFcenter(obj, indicesOfConesPooledByTheTargetRFcenter, spatialSupportDegs);

    % Compute the target visual RF map
    paramsVector(1) = 1;
    paramsVector(2) = targetVisualRFDoGparams.surroundToCenterRcRatio;
    paramsVector(3) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
    targetVisualRFmap = RetinaToVisualFieldTransformer.differenceOfArbitraryCenterAndGaussianSurroundRF(...
       visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs, paramsVector, spatialSupportDegs);
   
    % OR
    paramsVector(1) = 1;
    paramsVector(2) = visualRFcenterCharacteristicRadiusDegs;
    paramsVector(3) = targetVisualRFDoGparams.surroundToCenterRcRatio;
    paramsVector(4) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
  
    targetVisualRFmapIdeal = RetinaToVisualFieldTransformer.differenceOfGaussianCenterAndGaussianSurroundRF(...
       paramsVector, spatialSupportDegs);
   

    figure(11); clf;
    subplot(2,3,1);
    imagesc(spatialSupportDegs(:,1)*60, spatialSupportDegs(:,2)*60, visualRFcenterConeMap/max(visualRFcenterConeMap(:)));
    set(gca, 'CLim', [-1 1])
    title('cone-input based RF center')
    axis 'image'

    subplot(2,3,2)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), targetVisualRFmap);
    set(gca, 'CLim', [-1 1])
    axis 'image'
    title('target RF')
    colormap(gray)

    subplot(2,3,3)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), targetVisualRFmapIdeal);
    set(gca, 'CLim', [-1 1])
    axis 'image'
    title('target RF (ideal)')
    colormap(gray)

    profileCenter = sum(visualRFcenterConeMap,1);
    profileRF = sum(targetVisualRFmap,1);
    profileRFideal = sum(targetVisualRFmapIdeal,1);
    maxProfile = max([max(profileCenter) max(profileRF) max(profileRFideal)]);

    subplot(2,3,4)
    plot(spatialSupportDegs(:,1), profileCenter/maxProfile, 'r-');
    set(gca, 'YLim', [-1 1]);
    
    subplot(2,3,5)
    plot(spatialSupportDegs(:,1),  profileRF/max(profileRF), 'r-'); hold on
    plot(spatialSupportDegs(:,1),  profileRFideal/max(profileRFideal), 'b-');
    set(gca, 'YLim', [-1 1]);

    subplot(2,3,6)
    plot(spatialSupportDegs(:,1),  profileRFideal/maxProfile, 'b-');
    set(gca, 'YLim', [-1 1]);


    retinalRFparams = [];
    weightsComputeFunctionHandle = [];
    modelConstants = [];
end

function [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs] = analyzeRFcenter(obj, ...
    indicesOfConesPooledByTheTargetRFcenter, spatialSupportDegs)
    
    % Compute the visual RF center map
    visualRFcenterConeMap = computeVisualRFcenterMapFromInputConesInArcMinSupport(obj, indicesOfConesPooledByTheTargetRFcenter, spatialSupportDegs);
    
    % Compute its characteristic radius by fitting a flat top gaussian
    % ellipsoid to it.
    if (numel(indicesOfConesPooledByTheTargetRFcenter) == 2)
       % We are using a circularly summetric PSF, so force the
       % orientation to match the orientation of the 2 cones
       cone1RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheTargetRFcenter(1),:);
       cone2RFpos = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheTargetRFcenter(2),:);
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

function visualRFcenterConeMap = computeVisualRFcenterMapFromInputConesInArcMinSupport(obj, indicesOfConesPooledByTheTargetRFcenter, spatialSupportDegs)

    % Compute the retinal RF center cone map
    theConeCharacteristicRadiiDegs = 0.204*sqrt(2.0)*obj.theConeMosaic.coneApertureDiametersDegs(indicesOfConesPooledByTheTargetRFcenter);
    theConePositionsDegs = obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheTargetRFcenter,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);
    retinalRFcenterConeMap = retinalRFcenterConeMapFromPooledConeInputs(theConeCharacteristicRadiiDegs, theConePositionsDegs, spatialSupportDegs);

    % Compute the visual RF center cone map via convolution of the retinalRFcenterConeMap with the PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.theCircularPSFData.data, 'same');

    % Unit amplitude
    visualRFcenterConeMap = visualRFcenterConeMap ./ max(visualRFcenterConeMap(:));
end

function RF2D = retinalRFcenterConeMapFromPooledConeInputs(coneRc, conePos, spatialSupport)
    
    [X,Y] = meshgrid(spatialSupport(:,1), spatialSupport(:,2));

    conesNumPooledByRFcenter = size(conePos,1);
    for iCone = 1:conesNumPooledByRFcenter
  
        theConeApertureRF = exp(-((X-conePos(iCone,1))/coneRc(iCone)).^2) .* ...
                            exp(-((Y-conePos(iCone,2))/coneRc(iCone)).^2);
        
        if (iCone == 1)
            RF2D = theConeApertureRF;
        else
            RF2D = RF2D + theConeApertureRF;
        end
    end
    RF2D = RF2D / max(RF2D(:));

end