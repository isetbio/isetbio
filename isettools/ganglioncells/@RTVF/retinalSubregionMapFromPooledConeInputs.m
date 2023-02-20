function retinalSubregionConeMap = retinalSubregionMapFromPooledConeInputs(obj, ...
    conePosDegs, coneApertureAreas, osLengthAttenuationFactors, coneWeights)

    % Generate XY spatial position grid
    spatialSupportDegsX = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
    spatialSupportDegsY = obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs;
    [X,Y] = meshgrid(spatialSupportDegsX, spatialSupportDegsY);
    XY = [X(:) Y(:)];

    % Initialize subregion to all 0.
    % The extent of the map is equal to obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs
    % If cones fall outside of this extent they are not included and a
    % warning is generated
    retinalSubregionConeMap = X*0;

    conesNumPooled = size(conePosDegs,1);
    if (conesNumPooled == 0)
        return;
    end

    % Compute insertion coords for the cone apertures
    insertionCoords = cell(1, conesNumPooled);
    kernelSize = size(obj.coneApertureBlurKernel,1);
    parfor iCone = 1:conesNumPooled
        % Compute aperture map insertion coordinates
        dd = sum((bsxfun(@minus, XY, conePosDegs(iCone,:))).^2,2);
        [~,idx] = min(dd(:));
        xo = X(idx); yo = Y(idx);
        [~,row] = min(abs(yo - spatialSupportDegsY));
        [~,col] = min(abs(xo - spatialSupportDegsX));
        m = (kernelSize-1)/2;
        rr = row + (-m:1:m);
        cc = col + (-m:1:m);
        insertionCoords{iCone} = [cc(:) rr(:)];
    end % iCone

    retinalSubregionConeMap = X * 0;
    conesNotIncluded = 0;

    apertureScalingBasedOnQuantalEfficacy  = false;
    if (strcmp(obj.stfComputeMethod, RTVF.modeledSTFcomputeMethod))
        apertureScalingBasedOnQuantalEfficacy = true;
    end

    for iCone = 1:conesNumPooled
        % Extract aperture map insertion coordinates
        cc = insertionCoords{iCone}(:,1);
        rr = insertionCoords{iCone}(:,2);
        if (min(rr) < 1) || (min(cc) < 1) || (max(rr)>size(retinalSubregionConeMap,1)) || (max(cc) > size(retinalSubregionConeMap,2))
            conesNotIncluded = conesNotIncluded + 1;
            continue;
        end

        % The effective cone aperture is the cone apertur blur kernel (at
        % the cone's blur zone) multiplied by the aperture area of the individual cone
        % and by the os-length attenuation factor for the individual cone
        if (apertureScalingBasedOnQuantalEfficacy)
            effectiveConeApertureRF = obj.coneApertureBlurKernel * coneApertureAreas(iCone) * osLengthAttenuationFactors(iCone);
        else
            effectiveConeApertureRF = obj.coneApertureBlurKernel;
        end

        % Accumulate the subregion cone map
        retinalSubregionConeMap(rr,cc) = retinalSubregionConeMap(rr,cc) + coneWeights(iCone) * effectiveConeApertureRF;
    end
 
    if (conesNotIncluded > 0)
        fprintf(2,'RTVF.retinalSubregionMapFromPooledConeInputs:: %d of the %d cones pooled by the continuous model were NOT included in the actual subregion map because they fell outside of the spatial support.\n', conesNotIncluded, conesNumPooled);
    end

end