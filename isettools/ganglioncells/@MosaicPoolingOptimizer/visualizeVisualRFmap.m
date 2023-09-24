function visualizeVisualRFmap(theVisualRFmapStruct, retinalRGCRFposDegs, theAxes, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);

    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;
    ff = p.Results.withFigureFormat;
    spatialSupportRangeArcMin = tickSeparationArcMin * 4;

    % Find the coords of the RF max
    [~,idx] = max(abs(theVisualRFmapStruct.theRFmap(:)));
    [row,col] = ind2sub(size(theVisualRFmapStruct.theRFmap),idx);

    xShiftDegs = mean(theVisualRFmapStruct.spatialSupportDegsX) - theVisualRFmapStruct.spatialSupportDegsX(col);
    yShiftDegs = mean(theVisualRFmapStruct.spatialSupportDegsY) - theVisualRFmapStruct.spatialSupportDegsY(row);
    pixelSizeDegs = theVisualRFmapStruct.spatialSupportDegsX(2)-theVisualRFmapStruct.spatialSupportDegsX(1);
    xShiftSamples = sign(xShiftDegs)*round(abs(xShiftDegs)/pixelSizeDegs);
    yShiftSamples = sign(yShiftDegs)*round(abs(yShiftDegs)/pixelSizeDegs);

    visualRGCRFposDegs(1) = theVisualRFmapStruct.spatialSupportDegsX(col);
    visualRGCRFposDegs(2) = theVisualRFmapStruct.spatialSupportDegsY(row);

    xyLimsArcMin = spatialSupportRangeArcMin/2*[-1 1];
    xyTicksArcMin = 0:(tickSeparationArcMin):60;
    xyTicksArcMin = [-fliplr(xyTicksArcMin(2:end)) xyTicksArcMin];
    xyTicks = 0:(tickSeparationArcMin/60):30;
    xyTicks = [-fliplr(xyTicks(2:end)) xyTicks];

    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;

    % The cell position
    XLims = visualRGCRFposDegs(1) + [spatialSupportDegs(1) spatialSupportDegs(end)];
    YLims = visualRGCRFposDegs(2) + [spatialSupportDegs(1) spatialSupportDegs(end)];


    surroundBoostedRF = theVisualRFmapStruct.theRFmap;
    surroundBoostedRF(surroundBoostedRF<0) = surroundBoostedRF(surroundBoostedRF<0)*10;

    theRFprofileX = squeeze(sum(theVisualRFmapStruct.theRFmap,1));
    theRFprofileY = squeeze(sum(theVisualRFmapStruct.theRFmap,2));
    maxProfile = max([max(abs(theRFprofileX(:))) max(abs(theRFprofileY(:)))]);

    theRFprofileX = circshift(theRFprofileX ,xShiftSamples);
    theRFprofileY = circshift(theRFprofileY ,yShiftSamples);


    spatialSupportArcMinX = theVisualRFmapStruct.spatialSupportDegsX*60;
    spatialSupportArcMinX = spatialSupportArcMinX - mean(spatialSupportArcMinX);

    spatialSupportArcMinY = theVisualRFmapStruct.spatialSupportDegsY*60;
    spatialSupportArcMinY = spatialSupportArcMinY - mean(spatialSupportArcMinY);


    % Render the RF map
    MSreadyPlot.render2DVisualRF(theAxes{1,1}, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        theVisualRFmapStruct.theRFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, [], false, ff);

    MSreadyPlot.renderRFprofile(theAxes{1,2}, spatialSupportArcMinX, theRFprofileX/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ...
        'space, x (arc min)', ff);

    MSreadyPlot.renderRFprofile(theAxes{1,3}, spatialSupportArcMinY, theRFprofileY/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, false, gridlessLineWeightingFuncions, ...
        'space, y (arc min)', ff);

    % Render the surround boosted RFmap
    MSreadyPlot.render2DVisualRF(theAxes{1,4}, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        surroundBoostedRF,  ...
        xyTicks, XLims, YLims, reverseXDir, 'surround boosted (x10)', true, ff);

    drawnow;


    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';


    % The visual RF map
    pdfFileName = sprintf('InterpolatedRGC_VisualRFmap');
    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);

    hFig = figure(100); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    MSreadyPlot.render2DVisualRF(thePLOSAxes{1,1}, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        theVisualRFmapStruct.theRFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, [], false, ff);

    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The visual RF map profile X
    pdfFileName = sprintf('InterpolatedRGC_VisualRFprofileX');
    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);

    hFig = figure(200); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    MSreadyPlot.renderRFprofile(thePLOSAxes{1,1}, spatialSupportArcMinX, theRFprofileX/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ...
        'space, x (arc min)', ff);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The visual RF map profile Y
    pdfFileName = sprintf('InterpolatedRGC_VisualRFprofileY');
    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);

    hFig = figure(300); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    MSreadyPlot.renderRFprofile(thePLOSAxes{1,1}, spatialSupportArcMinY, theRFprofileY/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, false, gridlessLineWeightingFuncions, ...
        'space, y (arc min)', ff);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end