function [retinalSubregionConeMap, retinalSubregionConeMapFlatTop] = retinalSubregionConeMapFromPooledConeInputs(...
    theConeMosaic, theConeIndices, theConeWeights, spatialSupportDegs, flatTopSaturationLevel)
    
    [X,Y] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    XY = [X(:) Y(:)];

    conePosDegs = theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    conesNumPooled = size(conePosDegs,1);
    if (conesNumPooled == 0)
        retinalSubregionConeMap = X*0;
        retinalSubregionConeMapFlatTop = X*0;
        return;
    end

    oiResDegs = spatialSupportDegs(2,1) - spatialSupportDegs(1,1);
    theConeApertureBlurKernel = computeConeApertureBlurKernel(theConeMosaic, theConeIndices, oiResDegs);

    threshold = 0.1*max(theConeApertureBlurKernel(:));
    theConeApertureBlurKernelFlatTop = theConeApertureBlurKernel;
    theConeApertureBlurKernelFlatTop(theConeApertureBlurKernelFlatTop>threshold) = threshold;

    insertionCoords = cell(1, conesNumPooled);
    parfor iCone = 1:conesNumPooled
        % Compute aperture map insertion coordinates
        dd = sum((bsxfun(@minus, XY, conePosDegs(iCone,:))).^2,2);
        [~,idx] = min(dd(:));
        xo = X(idx); yo = Y(idx);
        [~,row] = min(abs(yo - spatialSupportDegs(:,2)));
        [~,col] = min(abs(xo - spatialSupportDegs(:,1)));
        m = (size(theConeApertureBlurKernel,1)-1)/2;
        rr = row + (-m:1:m);
        cc = col + (-m:1:m);
        insertionCoords{iCone} = [cc(:) rr(:)];
    end

    retinalSubregionConeMap = X * 0;
    retinalSubregionConeMapFlatTop = X*0;
    conesNotIncluded = 0;

    for iCone = 1:conesNumPooled
        % Extract aperture map insertion coordinates
        cc = insertionCoords{iCone}(:,1);
        rr = insertionCoords{iCone}(:,2);
        if (min(rr) < 1) || (min(cc) < 1) || (max(rr)>size(retinalSubregionConeMap,1)) || (max(cc) > size(retinalSubregionConeMap,2))
            conesNotIncluded = conesNotIncluded + 1;
            continue;
        end

        % Accumulate the subregion cone map
        retinalSubregionConeMap(rr,cc) = retinalSubregionConeMap(rr,cc) + theConeWeights(iCone) * theConeApertureBlurKernel;
        retinalSubregionConeMapFlatTop(rr,cc) = retinalSubregionConeMapFlatTop(rr,cc) + theConeWeights(iCone) * theConeApertureBlurKernelFlatTop;
    end

    if (conesNotIncluded > 0)
        fprintf(2,'%d of the %d cones pooled by the continuous model were NOT included in the actual subregion map because they fell outside of the spatial support.\n', conesNotIncluded, conesNumPooled);
    end

    saturationLevel = max(retinalSubregionConeMapFlatTop(:)) * flatTopSaturationLevel;
    idx = find(retinalSubregionConeMapFlatTop > saturationLevel);
    retinalSubregionConeMapFlatTop(idx) = saturationLevel;
end

function theConeApertureBlurKernel = computeConeApertureBlurKernel(theConeMosaic, theConeIndices, oiResDegs)

    % Find the blurZone in which theConeIndices belong to
    targetZoneIndex = coneApertureBlurZone(theConeMosaic, theConeIndices);

    % Retrieve the blurApertureDiameter for the target zone
    blurApertureDiameterDegs = theConeMosaic.blurApertureDiameterDegsZones(targetZoneIndex);

    % Compute the aperture blur kernel
    lowOpticalImageResolutionWarning = true;
    theConeApertureBlurKernel = theConeMosaic.generateApertureKernel(...
        blurApertureDiameterDegs(1), oiResDegs, lowOpticalImageResolutionWarning);

end

function targetZoneIndex = coneApertureBlurZone(theConeMosaic, theConeIndices)
    % Find the blur aperture diameter zone for theConeIndices 
    coneApertureBlurZonesNum = numel(theConeMosaic.blurApertureDiameterMicronsZones);
    inputConesNumInZone = zeros(1,coneApertureBlurZonesNum);

    for zoneIndex = 1:coneApertureBlurZonesNum 
        % Determine which cones should receive this blur.
        coneIDsInZone = theConeMosaic.coneIndicesInZones{zoneIndex};
        inputConesNumInZone(zoneIndex) = numel(find(ismember(theConeIndices,coneIDsInZone) == 1));
    end

    [~, targetZoneIndex] = max(inputConesNumInZone);
end