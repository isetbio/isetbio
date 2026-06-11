function [subregionLineWeightingFunctions, subregionContourData] = renderSubregionConePoolingPlot(ax, theConeMosaic, ...
        rgcRFposDegs, coneIndices, coneWeights, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTicks', false, @islogical);
    p.addParameter('noYTicks', false, @islogical);
    p.addParameter('gridless', false, @islogical);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('flatTopSaturationLevel', 1, @isscalar);
    p.addParameter('alsoComputeSubregionContour', false, @islogical);
    p.addParameter('overlayedSubregionContour', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('spatialSupportRangeArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('xAxisTickAngleRotationDegs', 0, @isscalar);
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('resetAxes', true, @islogical);

    p.parse(varargin{:});
    
    spatialSupportRangeArcMin = p.Results.spatialSupportRangeArcMin;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    computeSubregionContour = p.Results.alsoComputeSubregionContour;
    overlayedSubregionContour = p.Results.overlayedSubregionContour;
    flatTopSaturationLevel = p.Results.flatTopSaturationLevel;
    gridless = p.Results.gridless;

    plotTitle = p.Results.plotTitle;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTicks = p.Results.noXTicks;
    noYTicks = p.Results.noYTicks;
    ff = p.Results.withFigureFormat;
    resetAxes = p.Results.resetAxes;

    xAxisTickAngleRotationDegs = p.Results.xAxisTickAngleRotationDegs;
    
    if (isempty(tickSeparationArcMin))
        tickSeparationArcMin = 3;
    end

    if (isempty(spatialSupportRangeArcMin))
        spatialSupportRangeArcMin = 10;
    end
    
    if (isempty(tickSeparationArcMin))
        tickSeparationArcMin = 5.0;
    end

    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;
    spatialSupportXYDegs(:,1) = rgcRFposDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = rgcRFposDegs(2) + spatialSupportDegs;
    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = rgcRFposDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = rgcRFposDegs(2) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];

    [retinalSubregionConeMap, retinalSubregionConeMapFlatTop] = retinalSubregionConeMapFromPooledConeInputs(...
        theConeMosaic, coneIndices, coneWeights, spatialSupportXYDegs, flatTopSaturationLevel);

    if (computeSubregionContour)
        spatialSupportSamples = 64;
        coneApertureSizeSpecifierForRGCRFplotting = 'spacing based';
        switch (coneApertureSizeSpecifierForRGCRFplotting)
            case 'spacing based'
                coneRFradii = 0.6*0.5*theConeMosaic.coneRFspacingsDegs;
            case 'characteristic radius based'
                coneRFradii = ...
                    theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                    theConeMosaic.coneApertureDiametersDegs;
            otherwise
                error('Unknown apertureSizeSpecifierForRGCRFplotting: ''%s''.', coneApertureSizeSpecifierForRGCRFplotting);
        end

        conePositions = theConeMosaic.coneRFpositionsDegs(coneIndices,:);
        xSupport = spatialSupportXYDegs(:,1);
        ySupport = spatialSupportXYDegs(:,2);
        subregionContourData = mRGCMosaic.subregionOutlineContourFromPooledCones(...
             conePositions, coneRFradii, coneWeights, ...
             xSupport, ySupport, spatialSupportSamples);
    else
        subregionContourData = [];
    end


    if (resetAxes)
        cla(ax, 'reset');
    end

    imagesc(ax, spatialSupportXYDegs(:,1), spatialSupportXYDegs(:,2), retinalSubregionConeMapFlatTop);
    hold(ax, 'on');
    

    if (~isempty(overlayedSubregionContour))
        S.Vertices = overlayedSubregionContour.vertices;
        S.Faces = overlayedSubregionContour.faces;
        S.FaceVertexCData = [0.5 0.5 0.5];
        S.FaceColor = 'flat';
        S.EdgeColor = [0 0 0];
        S.FaceAlpha = 0.0;
        S.LineWidth = 2.0;
        patch(S, 'Parent', ax);
    end

    for iInputCone = 1:numel(coneIndices)
        switch theConeMosaic.coneTypes(coneIndices(iInputCone))
            case cMosaic.LCONE_ID
                coneColor = [1 0 0];
            case cMosaic.MCONE_ID
                coneColor = [0 0.75 0];
            case cMosaic.SCONE_ID
                coneColor = [0 0 1];
        end
        plot(ax,theConeMosaic.coneRFpositionsDegs(coneIndices(iInputCone),1), theConeMosaic.coneRFpositionsDegs(coneIndices(iInputCone),2), '.', 'MarkerSize', 14, 'Color', coneColor);
    end
    
    % Return the subregion line weighting functions along X and Y
    subregionLineWeightingFunctions.x = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,1), ...
        'amplitude', sum(retinalSubregionConeMap,1));

    subregionLineWeightingFunctions.y = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,2), ...
        'amplitude', sum(retinalSubregionConeMap,2));

    axis(ax, 'image'); axis(ax, 'xy');

    xyTicks = -30:(tickSeparationArcMin/60):30;
    if (tickSeparationArcMin >= 6)
        xyTickLabels = sprintf('%2.1f\n', xyTicks);
    else
        xyTickLabels = sprintf('%2.2f\n', xyTicks);
    end

    set(ax, 'CLim', [0 max(retinalSubregionConeMapFlatTop(:))], ...
            'XLim', XLims, 'YLim', YLims, ...
            'XTick', xyTicks, 'YTick', xyTicks, ...
            'XTickLabel', xyTickLabels, ...
            'YTickLabel', xyTickLabels);

    if (noXTicks)
        set(ax, 'XTickLabel', {});
    end

    if (noYTicks)
        set(ax, 'YTickLabel', {});
    end

    if (isempty(ff))
        if (~noXLabel)
            xlabel(ax, 'space, x (deg)');
        end
    
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)');
        end

        if (~isempty(plotTitle))
            title(ax, plotTitle);
        end

    else
        % Font size
        set(ax, 'FontSize', ff.fontSize);
    
        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
        if (~noXLabel)
            xlabel(ax, 'space, x (deg)', 'FontAngle', ff.axisFontAngle);
        end
    
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)' ,'FontAngle', ff.axisFontAngle);
        end

        if (~isempty(plotTitle))
            title(ax, plotTitle, ...
                'Color', ff.titleColor, 'FontSize', ff.titleFontSize, ...
                'FontWeight', ff.titleFontWeight);
        end

    end

    if (~gridless)
        grid(ax, 'on');
    else
        grid(ax, 'off');
    end

    xtickangle(ax, xAxisTickAngleRotationDegs);
    colormap(ax, brewermap(1024, 'greys'));
    
end

function [retinalSubregionConeMap, retinalSubregionConeMapFlatTop] = retinalSubregionConeMapFromPooledConeInputs(...
    theConeMosaic, theConeIndices, theConeWeights, spatialSupportDegs, flatTopSaturationLevel)
    
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

    halfMax = 0.1*max(theConeApertureBlurKernel(:));
    theConeApertureBlurKernelFlatTop = theConeApertureBlurKernel;
    theConeApertureBlurKernelFlatTop(theConeApertureBlurKernelFlatTop>halfMax) = halfMax;

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

