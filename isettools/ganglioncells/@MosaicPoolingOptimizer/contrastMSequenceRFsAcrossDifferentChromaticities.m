function contrastMSequenceRFsAcrossDifferentChromaticities(...
    exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
    theComputeReadyMRGCmosaic, ...
    theOptimallyMappedAchromaticRFmapsFileName, ...
    theOptimallyMappedLconeIsolatingRFmapsFileName, ...
    theOptimallyMappedMconeIsolatingRFmapsFileName, ...
    rfPixelsAcross, opticsParams, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript, ...
        varargin)

    % Parse input
    p = inputParser;
    p.addParameter('examinedRGCindices', [], @(x)(isnumeric(x)||(ischar(x))));
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.addParameter('performSurroundAnalysisForConesExclusiveToTheSurround', true, @islogical);
    p.addParameter('generateVideoWithAllExaminedRGCs', false, @islogical);

    p.parse(varargin{:});
    examinedRGCindices = p.Results.examinedRGCindices;
    performSurroundAnalysisForConesExclusiveToTheSurround = p.Results.performSurroundAnalysisForConesExclusiveToTheSurround;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    generateVideoWithAllExaminedRGCs = p.Results.generateVideoWithAllExaminedRGCs;


    % Load the achromatic RF maps
    load(theOptimallyMappedAchromaticRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');
    theOptimallyMappedAchromaticVisualRFmaps = optimallyMappedVisualRFmaps;
    clear 'optimallyMappedVisualRFmaps';

    % Load the L-cone isolating RF maps
    load(theOptimallyMappedLconeIsolatingRFmapsFileName, 'optimallyMappedVisualRFmaps');
    theOptimallyMappedLconeIsolatingVisualRFmaps = optimallyMappedVisualRFmaps;
    clear 'optimallyMappedVisualRFmaps';

    load(theOptimallyMappedMconeIsolatingRFmapsFileName, 'optimallyMappedVisualRFmaps');
    theOptimallyMappedMconeIsolatingVisualRFmaps = optimallyMappedVisualRFmaps;
    clear 'optimallyMappedVisualRFmaps';


    % Cone contrasts and overall contrast for achromatic stimulus
    [achromaticStimulusConeContrasts, achromaticStimulusContrast] = ...
        MosaicPoolingOptimizer.contrastForChromaticity('achromatic');

    % Cone contrasts and overall contrast for L-coneIsolating stimulus
    [LconeIsolatingStimulusConeContrasts, LconeIsolatingStimulusContrast] = ...
        MosaicPoolingOptimizer.contrastForChromaticity('Lcone isolating');

    % Cone contrasts and overall contrast for M-coneIsolating stimulus
    [MconeIsolatingStimulusConeContrasts, MconeIsolatingStimulusContrast] = ...
        MosaicPoolingOptimizer.contrastForChromaticity('Mcone isolating');

    % Correction factors to account for low L/M cone isolating stimulus contrast
    % compared to the achromatic sitmulus contrast
    LconeStimulusContrastCorrection =  achromaticStimulusContrast * achromaticStimulusConeContrasts(1) / (LconeIsolatingStimulusContrast * LconeIsolatingStimulusConeContrasts(1));
    MconeStimulusContrastCorrection =  achromaticStimulusContrast * achromaticStimulusConeContrasts(2) / (MconeIsolatingStimulusContrast * MconeIsolatingStimulusConeContrasts(2));


    if (isempty(examinedRGCindices))
        examinedRGCindices = 1:theComputeReadyMRGCmosaic.rgcsNum;
    end

    if (opticsParams.examinedSubjectRankOrder == 0)
       opticsString = sprintf('diffr-limited, pupil_%2.0fmm, defocus_%1.2fD', opticsParams.pupilDiameterMM, opticsParams.refractiveErrorDiopters);
    else
       opticsString = sprintf('physio-optics, defocus_%1.2fD', opticsParams.refractiveErrorDiopters);
    end

    analyzedLconeCenterRGCs = 0;
    analyzedMconeCenterRGCs = 0;

    for iRGC = 1:numel(examinedRGCindices)

        % Get theRGCindex
        theRGCindex = examinedRGCindices(iRGC);
        [~, ~, theCenterMajorityConeType] = theComputeReadyMRGCmosaic.centerConeTypeWeights(theRGCindex);

        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            analyzedLconeCenterRGCs = analyzedLconeCenterRGCs +1;
        end

        if (theCenterMajorityConeType == cMosaic.MCONE_ID)
            analyzedMconeCenterRGCs = analyzedMconeCenterRGCs +1;
        end
    end


    theLconeCenterData = cell(1, analyzedLconeCenterRGCs);
    theMconeCenterData = cell(1, analyzedMconeCenterRGCs);

    currentLconeCenterRGCindex = 0;
    currentMconeCenterRGCindex = 0;

    for iRGC = 1:numel(examinedRGCindices)

        % Get theRGCindex
        theRGCindex = examinedRGCindices(iRGC);
      
        % Analyze cone weights to center and surround
        [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix] = MosaicPoolingOptimizer.analyzeCenterSurroundConeMix(...
            theComputeReadyMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround);


        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            currentLconeCenterRGCindex = currentLconeCenterRGCindex + 1;
            theLconeCenterData{currentLconeCenterRGCindex} = struct(...
                'theRGCindex', theRGCindex, ...
                'netCenterLconeWeight', netCenterLconeWeight, ...
                'netCenterMconeWeight', netCenterMconeWeight, ...
                'netSurroundLconeWeight', netSurroundLconeWeight, ...
                'netSurroundMconeWeight', netSurroundMconeWeight, ...
                'surroundConeMix',  surroundConeMix, ...
                'achromaticRFmapDataStruct', theOptimallyMappedAchromaticVisualRFmaps{iRGC}, ...
                'LconeIsolatingRFmapDataStruct', theOptimallyMappedLconeIsolatingVisualRFmaps{iRGC}, ...
                'MconeIsolatingRFmapDataStruct', theOptimallyMappedMconeIsolatingVisualRFmaps{iRGC} ...
                );
        end

        if (theCenterMajorityConeType == cMosaic.MCONE_ID)
            currentMconeCenterRGCindex = currentMconeCenterRGCindex + 1;
            theMconeCenterData{currentMconeCenterRGCindex} = struct(...
                'theRGCindex', theRGCindex, ...
                'netCenterLconeWeight', netCenterLconeWeight, ...
                'netCenterMconeWeight', netCenterMconeWeight, ...
                'netSurroundLconeWeight', netSurroundLconeWeight, ...
                'netSurroundMconeWeight', netSurroundMconeWeight, ...
                'surroundConeMix',  surroundConeMix, ...
                'achromaticRFmapDataStruct', theOptimallyMappedAchromaticVisualRFmaps{iRGC}, ...
                'LconeIsolatingRFmapDataStruct', theOptimallyMappedLconeIsolatingVisualRFmaps{iRGC}, ...
                'MconeIsolatingRFmapDataStruct', theOptimallyMappedMconeIsolatingVisualRFmaps{iRGC} ...
                );
        end

    end

    plotSelectedRGCData(theComputeReadyMRGCmosaic, ...
        theLconeCenterData, theMconeCenterData, ...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        LconeStimulusContrastCorrection, MconeStimulusContrastCorrection, ...
        rfPixelsAcross, opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

end

function plotSelectedRGCData(theComputeReadyMRGCmosaic, ...
        theLconeCenterData, theMconeCenterData, ...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        LconeStimulusContrastCorrection, MconeStimulusContrastCorrection, ...
        rfPixelsAcross, opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    lConeCenterRGCsNum = numel(theLconeCenterData);
    mConeCenterRGCsNum = numel(theMconeCenterData);

    theLconeCenterRFpositions = zeros(lConeCenterRGCsNum ,2);
    theMconeCenterRFpositions = zeros(mConeCenterRGCsNum ,2);

    for iLconeRGC = 1:lConeCenterRGCsNum
        d = theLconeCenterData{iLconeRGC};
        theLconeCenterRFpositions(iLconeRGC,:) = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(d.theRGCindex,:);
    end

    for iMconeRGC = 1:mConeCenterRGCsNum
        d = theMconeCenterData{iMconeRGC};
        theMconeCenterRFpositions(iMconeRGC,:) = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(d.theRGCindex,:);
    end

    distances = sum((bsxfun(@minus, theLconeCenterRFpositions, exampleLconeCenterRGCposition)).^2,2);
    [~,idx] = min(distances(:));
    if (~isempty(idx))
        exampleLconeCenterData = theLconeCenterData{idx};
    else
        exampleLconeCenterData = [];
    end

    distances = sum((bsxfun(@minus, theMconeCenterRFpositions, exampleMconeCenterRGCposition)).^2,2);
    [~, idx] = min(distances(:));
    if (~isempty(idx))
        exampleMconeCenterData = theMconeCenterData{idx};
    else
        exampleMconeCenterData = [];
    end

    
    
    renderMSequenceRFMaps(20, theComputeReadyMRGCmosaic, exampleLconeCenterData, 'L', ...
        LconeStimulusContrastCorrection, MconeStimulusContrastCorrection, ...
        rfPixelsAcross, opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);

    renderMSequenceRFMaps(21, theComputeReadyMRGCmosaic, exampleMconeCenterData, 'M', ...
        LconeStimulusContrastCorrection, MconeStimulusContrastCorrection, ...
        rfPixelsAcross, opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);

end

function renderMSequenceRFMaps(figNo, theComputeReadyMRGCmosaic, exampleConeCenterData, coneType, ...
    LconeStimulusContrastCorrection, MconeStimulusContrastCorrection, ...
    rfPixelsAcross, opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    % Apply correction factors to account for low L/M cone isolating stimulus contrast
    % compared to the achromatic stimulus contrast
    exampleConeCenterData.LconeIsolatingRFmapDataStruct.theRFmap = ...
        exampleConeCenterData.LconeIsolatingRFmapDataStruct.theRFmap * LconeStimulusContrastCorrection;
    exampleConeCenterData.MconeIsolatingRFmapDataStruct.theRFmap = ...
        exampleConeCenterData.MconeIsolatingRFmapDataStruct.theRFmap * MconeStimulusContrastCorrection;


    % Smoothed maps
    [exampleConeCenterData.achromaticRFmapDataStruct.theSmoothedRFmap, smoothingKernel, rfPixelSizeSamples] = MosaicPoolingOptimizer.applyReidShapleySmoothingToRFmap(...
       exampleConeCenterData.achromaticRFmapDataStruct.spatialSupportDegsX, exampleConeCenterData.achromaticRFmapDataStruct.theRFmap, rfPixelsAcross);

    exampleConeCenterData.LconeIsolatingRFmapDataStruct.theSmoothedRFmap = MosaicPoolingOptimizer.applyReidShapleySmoothingToRFmap(...
       exampleConeCenterData.LconeIsolatingRFmapDataStruct.spatialSupportDegsX, exampleConeCenterData.LconeIsolatingRFmapDataStruct.theRFmap, rfPixelsAcross);

    exampleConeCenterData.MconeIsolatingRFmapDataStruct.theSmoothedRFmap = MosaicPoolingOptimizer.applyReidShapleySmoothingToRFmap(...
       exampleConeCenterData.MconeIsolatingRFmapDataStruct.spatialSupportDegsX, exampleConeCenterData.MconeIsolatingRFmapDataStruct.theRFmap, rfPixelsAcross);

    % Fit ellipse to the achromatic RF
    [achromaticRFcenterContourData, achromaticRFcentroid] = ellipseContourFromSubregionRFmap(exampleConeCenterData.achromaticRFmapDataStruct);
    
    

    % Max range
    m = max([max(abs(exampleConeCenterData.LconeIsolatingRFmapDataStruct.theRFmap(:))) max(abs(exampleConeCenterData.MconeIsolatingRFmapDataStruct.theRFmap(:)))]);
    RFmapRange = 0.5 * m * [-1 1];


    [cLUToriginal, cLUTreversed]= MosaicPoolingOptimizer.generateReidShapleyRFmapLUT();
    cLUT = cLUToriginal;
    cLUT = brewermap(1025, '*RdBu');

    visualizedRangeDegs = 0.2;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 800 800]);

    ax = subplot(3,3,1);
    renderRFmap(ax, exampleConeCenterData.achromaticRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData, achromaticRFcentroid, visualizedRangeDegs, cLUT, 'raw', 'achromatic RF');

    ax = subplot(3,3,2);
    renderRFmap(ax, exampleConeCenterData.LconeIsolatingRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData, achromaticRFcentroid, visualizedRangeDegs, cLUT, 'raw', 'L-cone RF');

    ax = subplot(3,3,3);
    renderRFmap(ax, exampleConeCenterData.MconeIsolatingRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData, achromaticRFcentroid, visualizedRangeDegs, cLUT, 'raw', 'M-cone RF');
    
    ax = subplot(3,3,4);
    renderRFmap(ax, exampleConeCenterData.achromaticRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData, achromaticRFcentroid, visualizedRangeDegs, cLUT, 'smoothed', 'achromatic RF');

    ax = subplot(3,3,5);
    renderRFmap(ax, exampleConeCenterData.LconeIsolatingRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData, achromaticRFcentroid, visualizedRangeDegs, cLUT, 'smoothed', 'L-cone RF');

    ax = subplot(3,3,6);
    renderRFmap(ax, exampleConeCenterData.MconeIsolatingRFmapDataStruct, RFmapRange, ...
        achromaticRFcenterContourData,  achromaticRFcentroid, visualizedRangeDegs, cLUT, 'smoothed', 'M-cone RF');
    
    ax = subplot(3,3,7);
    YLims = [];
    YLims = renderRadiallySymmetricRFprofile(ax, exampleConeCenterData.achromaticRFmapDataStruct, YLims, ...
        achromaticRFcentroid, achromaticRFcentroid, visualizedRangeDegs);

    ax = subplot(3,3,8);
    renderRadiallySymmetricRFprofile(ax, exampleConeCenterData.LconeIsolatingRFmapDataStruct, YLims, ...
        [], achromaticRFcentroid, visualizedRangeDegs);

    ax = subplot(3,3,9);
    renderRadiallySymmetricRFprofile(ax, exampleConeCenterData.MconeIsolatingRFmapDataStruct, YLims, ...
        [], achromaticRFcentroid, visualizedRangeDegs);

    drawnow;
end


function YLims = renderRadiallySymmetricRFprofile(ax, RFmapDataStruct, YLims, theRFmapCentroid, theRFmapReferenceCentroid, visualizedRangeDegs)

    if (isempty(theRFmapCentroid))
        [~, theRFmapCentroid] = ellipseContourFromSubregionRFmap(RFmapDataStruct);
    end

    [~,theRFcenterCol] = min(abs(RFmapDataStruct.spatialSupportDegsX-theRFmapCentroid(1)));
    [~,theRFcenterRow] = min(abs(RFmapDataStruct.spatialSupportDegsY-theRFmapCentroid(2)));

    radiallySymmetricRFprofileX = circularlyAveragedRFmap(RFmapDataStruct.theSmoothedRFmap, theRFcenterCol, theRFcenterRow);

    faceColor = [0.8 0.8 0.8];
    edgeColor = [0.2 0.2 0.2];
    faceAlpha = 0.8;
    lineWidth = 1.0;
    shadedAreaPlot(ax,RFmapDataStruct.spatialSupportDegsX, radiallySymmetricRFprofileX, radiallySymmetricRFprofileX*0, ...
        faceColor, edgeColor, faceAlpha, lineWidth)
    
    hold(ax, 'on');
    if (isempty(YLims))
        YLims(2) = 1.3*max(abs(radiallySymmetricRFprofileX(:)));
        YLims(1) =-YLims(2)*0.2;
    end
    YTicks = linspace(0, YLims(2), 5);
    YTicks = [-fliplr(YTicks) YTicks(2:end)];
    XTicks = theRFmapReferenceCentroid(1) + (-2:0.05:2);

    plot(ax, theRFmapReferenceCentroid(1)*[1 1], YLims, 'k-', 'LineWidth', 1.0);
    plot(ax, visualizedRangeDegs*0.5*[-1 1] + theRFmapReferenceCentroid(1), [0 0], 'k-', 'LineWidth', 1.0);
    box(ax, 'off');
    grid(ax, 'on');
    axis(ax, 'square');
    set(ax, 'XLim', visualizedRangeDegs*0.5*[-1 1] + theRFmapReferenceCentroid(1), 'YLim', YLims);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'XTick', XTicks, 'YTick', YTicks, 'XTickLabel', sprintf('%2.2f\n', XTicks), 'YTickLabel', {}); 
    set(ax, 'FontSize', 16);
end

function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end

function [theRadiallySymmetricRFprofileX, theRadiallySymmetricRFprofileY] = circularlyAveragedRFmap(theRF, theRFcenterCol, theRFcenterRow)
    % Compute unit circle radial distances
    N = size(theRF,1);
    xx = (1:N)-theRFcenterCol;
    yy = (1:N)-theRFcenterRow;
    [X,Y] = meshgrid(xx,yy);
    r = sqrt(X.^2+Y.^2);

    % Compute bins of data matrix for averaging by accumarray
    binLocations = round(r)+1;

    % reshape matrices to vectors for accumarray operation
    binLocations = reshape(binLocations, N*N, 1);
    theRF = reshape(theRF, N*N, 1);

    % Compute radial average (mean)
    theRadialRFprofile = accumarray(binLocations,theRF,[],@mean);
  
    nn = 1:(N+1-theRFcenterCol);
    theRadiallySymmetricRFprofileX = zeros(1,N);
    for i = 1:numel(nn)
        if (theRFcenterCol-1+nn(i) <= N) && (theRFcenterCol-1+nn(i) >= 1)
            theRadiallySymmetricRFprofileX(theRFcenterCol-1+nn(i)) = theRadialRFprofile(nn(i));
        end
    end
    for i = 1:numel(nn)
        if (theRFcenterCol+1-nn(i) <= N) && (theRFcenterCol+1-nn(i) >= 1)
            theRadiallySymmetricRFprofileX(theRFcenterCol+1-nn(i)) = theRadialRFprofile(nn(i));
        end
    end


    nn = 1:(N+1-theRFcenterRow);
    theRadiallySymmetricRFprofileY = zeros(1,N);
    for i = 1:numel(nn)
        if (theRFcenterRow-1+nn(i) <= N) && (theRFcenterRow-1+nn(i) >= 1)
            theRadiallySymmetricRFprofileY(theRFcenterRow-1+nn(i)) = theRadialRFprofile(nn(i));
        end
    end
    for i = 1:numel(nn)
        if (theRFcenterRow+1-nn(i) <= N) && (theRFcenterRow+1-nn(i) >= 1)
            theRadiallySymmetricRFprofileY(theRFcenterRow+1-nn(i)) = theRadialRFprofile(nn(i));
        end
    end

end




function renderRFmap(ax, RFmapDataStruct, RFmapRange, contourData, RFcentroid, visualizedRangeDegs, cLUT, whichMap, plotTitle)

    switch (whichMap)
        case 'raw'
            imagesc(ax, RFmapDataStruct.spatialSupportDegsX, ...
                        RFmapDataStruct.spatialSupportDegsY, ...
                        RFmapDataStruct.theRFmap);
            
        case 'smoothed'
            imagesc(ax, RFmapDataStruct.spatialSupportDegsX, ...
                        RFmapDataStruct.spatialSupportDegsY, ...
                        RFmapDataStruct.theSmoothedRFmap);
    end
    hold(ax, 'on');
    if (~isempty(contourData))
        S = contourData;
        S.FaceVertexCData = [0.85 0.85 0.85];
            S.FaceColor = 'flat';
            S.EdgeColor = [0 0 0];
            S.FaceAlpha = 0.0;
            S.LineWidth = 1.0;
            patch(S, 'Parent', ax)

    end

    XTicks = RFcentroid(1) + (-2:0.05:2);
    YTicks = RFcentroid(2) + (-2:0.05:2);

    plot(ax, RFcentroid(1)*[1 1], [RFmapDataStruct.spatialSupportDegsY(1) RFmapDataStruct.spatialSupportDegsY(end)], 'k-', 'LineWidth', 1.0);
    plot(ax,  [RFmapDataStruct.spatialSupportDegsX(1) RFmapDataStruct.spatialSupportDegsX(end)], RFcentroid(2)*[1 1], 'k-', 'LineWidth', 1.0);
    hold(ax, 'off');
    colormap(ax, cLUT);
    midPoint = (size(cLUT,1)-1)/2+1;
    %colorbar
    set(ax, 'Color', cLUT(midPoint,:));
    set(ax, 'CLim', RFmapRange);
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'XLim', visualizedRangeDegs*0.5*[-1 1] + RFcentroid(1));
    set(ax, 'YLim', visualizedRangeDegs*0.5*[-1 1] + RFcentroid(2));
    set(ax, 'XTick', XTicks, 'YTick', YTicks, 'XTickLabel', sprintf('%2.2f\n', XTicks), 'YTickLabel', sprintf('%2.2f\n', YTicks));
    set(ax, 'FontSize', 16);
    title(ax, plotTitle);
end


function [contourData, centroidPos] = ellipseContourFromSubregionRFmap(RFmapDataStruct)

    xSupport = RFmapDataStruct.spatialSupportDegsX;
    ySupport = RFmapDataStruct.spatialSupportDegsY;
   

    zLevel = 0.1;
    centerSubregionContourSamples = 30;

    RF = RFmapDataStruct.theSmoothedRFmap;
    maxRF = max(RF(:));
    minRF = min(RF(:));
    if (abs(minRF) > maxRF)
        RF = -RF;
    end

    % Binarize
    RF = RF / max(RF(:));
    RF(RF<zLevel) = 0.0;
    RF(RF>0) = 1.0;
    BW = imbinarize(RF);

    % Extract the maximum area
    BW = imclearborder(BW);
    BW = bwareafilt(BW,1);

    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    if (isempty(s))
        fprintf(2, 'Could not fit an ellipse to subregion map. Returning empty contour.\n')
        contourData = [];
        return;
    end
    
    % Calculate the ellipse line
    theta = linspace(0, 2*pi, centerSubregionContourSamples);
    col = (s.MajorAxisLength/2)*cos(theta);
    row = (s.MinorAxisLength/2)*sin(theta);
    M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];

    x = D(1,:);
    y = D(2,:);
    x = (x-1)/(numel(xSupport)-1) * (xSupport(end)-xSupport(1)) + xSupport(1); 
    y = (y-1)/(numel(ySupport)-1) * (ySupport(end)-ySupport(1)) + ySupport(1); 

    centroidPos(1) = mean(x(:));
    centroidPos(2) = mean(y(:));

    v = [x(:) y(:)];
    f = 1:numel(x);
    s = struct('faces', f, 'vertices', v);
    contourData = s;
    
end


