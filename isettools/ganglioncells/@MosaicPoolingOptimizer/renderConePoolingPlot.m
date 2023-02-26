function lineWeightingFunctions = renderConePoolingPlot(ax, theConeMosaic, ...
        rgcRFposDegs, coneIndices, coneWeights)

    spatialSupportDegs = (-20:0.1:20)/60;
    spatialSupportXYDegs(:,1) = rgcRFposDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = rgcRFposDegs(2) + spatialSupportDegs;
    XLims = rgcRFposDegs(1) + [spatialSupportDegs(1) spatialSupportDegs(end)];
    YLims = rgcRFposDegs(2) + [spatialSupportDegs(1) spatialSupportDegs(end)];

    retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(...
        theConeMosaic, coneIndices, coneWeights, spatialSupportXYDegs);

    imagesc(ax, spatialSupportXYDegs(:,1), spatialSupportXYDegs(:,2), retinalSubregionConeMap);
    hold(ax, 'on');
    for iInputCone = 1:numel(coneIndices)
        switch theConeMosaic.coneTypes(coneIndices(iInputCone))
            case cMosaic.LCONE_ID
                coneColor = [1 0 0];
            case cMosaic.MCONE_ID
                coneColor = [0 1 0];
            case cMosaic.SCONE_ID
                coneColor = [0 0 1];
        end
        plot(ax,theConeMosaic.coneRFpositionsDegs(coneIndices(iInputCone),1), theConeMosaic.coneRFpositionsDegs(coneIndices(iInputCone),2), '.', 'MarkerSize', 8, 'Color', coneColor);
    end
    

    lineWeightingFunctions.x = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,1), ...
        'amplitude', sum(retinalSubregionConeMap,1));

    lineWeightingFunctions.y = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,2), ...
        'amplitude', sum(retinalSubregionConeMap,2));
    axis(ax, 'image');
    axis(ax, 'xy');
    set(ax, 'CLim', [0 0.3*max(retinalSubregionConeMap(:))], ...
            'XLim', XLims, 'YLim', YLims);

    colormap(ax, brewermap(1024, 'greys'));
end

function retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(...
    theConeMosaic, theConeIndices, theConeWeights, spatialSupportDegs)
    
    [X,Y] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    XY = [X(:) Y(:)];

    conePosDegs = theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    conesNumPooled = size(conePosDegs,1);
    if (conesNumPooled == 0)
        retinalSubregionConeMap = X*0;
        return;
    end

    oiResDegs = spatialSupportDegs(2,1) - spatialSupportDegs(1,1);
    theConeApertureBlurKernel = computeConeApertureBlurKernel(theConeMosaic, theConeIndices, oiResDegs);

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
    end

    if (conesNotIncluded > 0)
        fprintf(2,'%d of the %d cones pooled by the continuous model were NOT included in the actual subregion map because they fell outside of the spatial support.\n', conesNotIncluded, conesNumPooled);
    end

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

