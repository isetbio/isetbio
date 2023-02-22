function [theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
    visualRFfromRetinalConePooling(obj, modelConstants, retinalPoolingParams)

    % Compute the indices and retinal weights of cones pooled by the center and the
    % surround mechanism
    pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(modelConstants, retinalPoolingParams);

    
    if (strcmp(obj.stfComputeMethod, RTVF.modeledSTFcomputeMethod))
        apertureScalingBasedOnQuantalEfficacy = true;
    else
        apertureScalingBasedOnQuantalEfficacy = false;
    end

    % Retrieve the center cone indices, their aperture areas and OS
    % length attenutation factors
    theCenterConeIndices = pooledConeIndicesAndWeights.centerConeIndices;
    if (apertureScalingBasedOnQuantalEfficacy)
        theConeApertureAreas = modelConstants.theConeMosaic.computeApertureAreasMetersSquared(theCenterConeIndices) * 1e16;
        theEffectiveOSlengthAttenuationFactors = modelConstants.theConeMosaic.computeEffectiveOSlengthAttenuationFactors(theCenterConeIndices);
    else
        theConeApertureAreas = [];
        theEffectiveOSlengthAttenuationFactors = [];
    end

    % Re-center the center cone positions
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(theCenterConeIndices ,:);
    rfCenterPositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);

    % Compute the retinal cone input map for the RF center
    theRetinalRFcenterConeMap = obj.retinalSubregionMapFromPooledConeInputs(...
        theConePositionsDegs, ...
        theConeApertureAreas, ...
        theEffectiveOSlengthAttenuationFactors, ...
        pooledConeIndicesAndWeights.centerConeWeights);

    % Retrieve the surround cone indices, their aperture areas and OS
    % length attenutation factors
    theSurroundConeIndices = pooledConeIndicesAndWeights.surroundConeIndices;
    
    if (apertureScalingBasedOnQuantalEfficacy)
        theConeApertureAreas = modelConstants.theConeMosaic.computeApertureAreasMetersSquared(theSurroundConeIndices) * 1e16;
        theEffectiveOSlengthAttenuationFactors = modelConstants.theConeMosaic.computeEffectiveOSlengthAttenuationFactors(theSurroundConeIndices);
    else
        theConeApertureAreas = [];
        theEffectiveOSlengthAttenuationFactors = [];
    end

    % Re-center the surround cone positions
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(theSurroundConeIndices,:);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);

    % Compute the retinal cone input map for the RF surround
    theRetinalRFsurroundConeMap  = obj.retinalSubregionMapFromPooledConeInputs(...
        theConePositionsDegs, ...
        theConeApertureAreas, ...
        theEffectiveOSlengthAttenuationFactors, ...
        pooledConeIndicesAndWeights.surroundConeWeights);

    if (apertureScalingBasedOnQuantalEfficacy)
        % Convolve the retinal RF center with the center cone PSF to get the corresponding visual RF center
        theVisualRFcenter = conv2(theRetinalRFcenterConeMap, modelConstants.theRFCenterConeMajorityPSF, 'same');
    
        % Convolve the retinal RF surround with the surround L+M PSF to get the corresponding visual RF surround
        theVisualRFsurround = conv2(theRetinalRFsurroundConeMap, modelConstants.theSurroundLconePlusMconePSF, 'same');
    
        % The composite visual RF cone map
        theVisualRF = theVisualRFcenter - theVisualRFsurround;
    else
        theVisualRF = [];
    end
end

