function retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(...
    theConeMosaic, conePosDegs, coneWeights, spatialSupportDegs)
    

    [X,Y] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    XY = [X(:) Y(:)];

    conesNumPooled = size(conePosDegs,1);
    if (conesNumPooled == 0)
        retinalSubregionConeMap = X*0;
        return;
    end

    dx = spatialSupportDegs(2,1) - spatialSupportDegs(1,1);
    dy = spatialSupportDegs(2,2) - spatialSupportDegs(1,2);

    lowOpticalImageResolutionWarning = true;
    oiResDegs = dx;
    zoneIndex = 1;
    if (numel(theConeMosaic.blurApertureDiameterDegsZones)>1)
        error('RTVFT.retinalSubregionConeMapFromPooledConeInputs(): the cone mosaic has more than 1 blur diameter zones. Dont know how to deal with this.')
    end

    theConeApertureRF = cell(1, conesNumPooled);
    insertionCoords = cell(1, conesNumPooled);
    parfor iCone = 1:conesNumPooled
       
        % Compute cone aperture map
        blurApertureDiameterDegs = theConeMosaic.blurApertureDiameterDegsZones(zoneIndex);
        theConeApertureRF{iCone} = theConeMosaic.generateApertureKernel(blurApertureDiameterDegs(1), oiResDegs, lowOpticalImageResolutionWarning);

        % Divide by relative efficacy because during actual computation  this relative cone efficacy is applied
        theConeApertureRF{iCone} = theConeApertureRF{iCone} ; 

        % Compute aperture map insertion coordinates
        dd = sum((bsxfun(@minus, XY, conePosDegs(iCone,:))).^2,2);
        [~,idx] = min(dd(:));
        xo = X(idx); yo = Y(idx);
        [~,row] = min(abs(yo - spatialSupportDegs(:,2)));
        [~,col] = min(abs(xo - spatialSupportDegs(:,1)));
        m = (size(theConeApertureRF{iCone})-1)/2;
        rr = row + (-m:1:m);
        cc = col + (-m:1:m);
        insertionCoords{iCone} = [cc(:) rr(:)];
    end

    retinalSubregionConeMap = X * 0;
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
        retinalSubregionConeMap(rr,cc) = retinalSubregionConeMap(rr,cc) + coneWeights(iCone) * theConeApertureRF{iCone};
    end

    if (conesNotIncluded > 0)
        fprintf(2,'%d of the %d cones pooled by the continuous model were NOT included in the actual subregion map because they fell outside of the spatial support.\n', conesNotIncluded, conesNumPooled);
        pause(1);
    end

end