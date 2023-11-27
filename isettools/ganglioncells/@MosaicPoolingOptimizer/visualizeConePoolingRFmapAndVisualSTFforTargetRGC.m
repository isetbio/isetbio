function visualizeConePoolingRFmapAndVisualSTFforTargetRGC(...
    computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename, pdfFileName, ...
    targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;

    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');
    load(mRGCMosaicSTFresponsesFilename, 'spatialFrequenciesTested', 'theMRGCMosaicOptimalSTFs');

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x4 RF poster');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    % Plot the retinal cone pooling RF map
    theVisualizedRGCindex = theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapNearPosition(...
        targetRGCposition, targetCenterConesNum, ...
        targetCenterConeMajorityType, ...
        'theAxes', theAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
        'withFigureFormat', ff, ...
        'reverseXDir', reverseXDir, ...
        'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions);

    theVisualizedRGCvisualSTFdata = theMRGCMosaicOptimalSTFs{theVisualizedRGCindex};
    idx = find(strcmp(theVisualizedRGCvisualSTFdata.DoGfitParams.names, 'kS/kC'));
    KsKcRatio = theVisualizedRGCvisualSTFdata.DoGfitParams.finalValues(idx);
    idx = find(strcmp(theVisualizedRGCvisualSTFdata.DoGfitParams.names, 'RsToRc'));
    RsRcRatio = theVisualizedRGCvisualSTFdata.DoGfitParams.finalValues(idx);
    SCintSensRatio = KsKcRatio * (RsRcRatio)^2;


    % Now plot the visual STF
    plotTitle = sprintf('R_{srnd}/R_{cntr} = %2.2f, S_{srnd}/S_{cntr} = %2.2f', RsRcRatio, SCintSensRatio);
    MSreadyPlot.renderSTF(theAxes{1,4}, ...
         spatialFrequenciesTested, theVisualizedRGCvisualSTFdata.measured, ...
         theVisualizedRGCvisualSTFdata.DoGfit.sfHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.compositeSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.centerSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.surroundSTFHiRes, ...
         plotTitle, ...
         [], ff, ...
         'noYLabel', ~true, ...
         'noYTickLabel', ~true, ...
         'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);

    % Export to PDF
    if (isnan(targetCenterConeMajorityType))
        pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM', ...
                    targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
    else
        if (isempty(targetCenterConeMajorityType))
            [~, theCenterConeTypeNum, targetCenterConeMajorityType, ~] = theComputeReadyMRGCmosaic.centerConeTypeWeights(theVisualizedRGCindex);
            targetCenterConesNum = sum(theCenterConeTypeNum);
        end
        if (isempty(targetRGCposition))
            targetRGCposition = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(theVisualizedRGCindex,:);
        end

        if (isnan(targetCenterConeMajorityType))
            pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);
        else
            switch (targetCenterConeMajorityType)
               case cMosaic.LCONE_ID
                    pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_LconeDominated_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);
               case cMosaic.MCONE_ID
                    pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_MconeDominated_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);                    
            end
        end
    end

    pdfFileName = strrep(pdfFileName, '.pdf', pdfPostFix);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300); 



    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';


    % The visual STF
    pdfFileName = sprintf('InterpolatedRGC_VisualSTF');

    hFigSTF = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

    MSreadyPlot.renderSTF(ax, ...
         spatialFrequenciesTested, theVisualizedRGCvisualSTFdata.measured, ...
         theVisualizedRGCvisualSTFdata.DoGfit.sfHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.compositeSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.centerSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.surroundSTFHiRes, ...
         '', ...
         [], ff, ...
         'noYLabel', ~true, ...
         'noYTickLabel', ~true, ...
         'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);
    text(0.11, 0.14, sprintf('\\rho = %2.2f', RsRcRatio), 'FontSize', ff.fontSize-4);
    text(0.11, 0.07, sprintf('\\sigma = %2.2f', SCintSensRatio), 'FontSize', ff.fontSize-4);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigSTF, 300);


    % The RF center cone pooling weights
    pdfFileNameSurroundWeightsNum = sprintf('InterpolatedRGC_RFsurroundWeightsMap');
    hFigRetinalConePoolingMap = figure(2); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigRetinalConePoolingMap,ff);
    set(hFigRetinalConePoolingMap, 'Color', [1 1 1]);
    axSurroundWeightsMap = theAxes{1,1};

    % The X-line weighting functions
    pdfFileNameXineWeightingFunctions  = sprintf('InterpolatedRGC_LineWeightingFunctionX');
    hFigLineWeightingFunctionsX = figure(3); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigLineWeightingFunctionsX,ff);
    set(hFigLineWeightingFunctionsX, 'Color', [1 1 1]);
    axLineWeightingFunctionsX = theAxes{1,1};

    % The Y-line weighting functions
    pdfFileNameYineWeightingFunctions = sprintf('InterpolatedRGC_LineWeightingFunctionY');
    hFigLineWeightingFunctionsY = figure(4); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigLineWeightingFunctionsY,ff);
    set(hFigLineWeightingFunctionsY, 'Color', [1 1 1]);
    axLineWeightingFunctionsY = theAxes{1,1};

    theAxes{1,1} = axSurroundWeightsMap;
    theAxes{1,2} = axLineWeightingFunctionsX;
    theAxes{1,3} = axLineWeightingFunctionsY;

    % Plot the retinal cone pooling RF map
    theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapOfRGCwithIndex(...
        theVisualizedRGCindex, ...
        'theAxes', theAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
        'withFigureFormat', ff, ...
        'reverseXDir', reverseXDir, ...
        'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions);


    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileNameSurroundWeightsNum);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigRetinalConePoolingMap, 300);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileNameXineWeightingFunctions);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigLineWeightingFunctionsX, 300);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileNameYineWeightingFunctions);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigLineWeightingFunctionsY, 300);

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);


end