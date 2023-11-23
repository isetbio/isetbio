function performContrastSTFsAcrossDifferentChromaticities(...
            mosaicParams, rawFiguresRoot, scaledFiguresRoot, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('performSurroundAnalysisForConesExclusiveToTheSurround', true, @islogical);
    p.addParameter('targetRangeForSurroundConeMix', [0.4 0.5], @(x)(isnumeric(x)&&(numel(x)==2)));
    p.parse(varargin{:});
    performSurroundAnalysisForConesExclusiveToTheSurround = p.Results.performSurroundAnalysisForConesExclusiveToTheSurround;
    targetRangeForSurroundConeMix = p.Results.targetRangeForSurroundConeMix;
    
    % Get PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
        'mosaicParams', mosaicParams);

    % Ask the user what optics were used for computing the compute-ready MRGC mosaic
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParamsForComputeReadyMosaic = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index was used for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(...
        mosaicParams, opticsParamsForComputeReadyMosaic);

    % Generate the filename of the compute-ready mRGCMosaic to use
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParamsForComputeReadyMosaic, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Now, ask the user what optics were used for computing the mRGC mosaic STF responses, so we can obtain the corresponding mRGCMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the mRGC mosaic STF responses \n');
    opticsParamsForMRGCSTFs = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate filename for the computed mRGCMosaicSTF responses
    [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForMRGCSTFs);

    % Achromatic STF responses filename
    coneFundamentalsOptimizedForStimPosition = false;
    [~, ~, mRGCMosaicAchromaticSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSTFresponsesFileName, 'STFresponses', ...
        'doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', struct(...
                    'stimulusChromaticity', 'achromatic',...
                    'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition));

    % L-cone isolating STF responses filename
    [~, ~, mRGCMosaicLconeIsolatingSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSTFresponsesFileName, 'STFresponses', ...
        'doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', struct(...
                    'stimulusChromaticity', 'Lcone isolating',...
                    'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition));

    % M-cone isolating STF responses filename
    [~, ~, mRGCMosaicMconeIsolatingSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSTFresponsesFileName, 'STFresponses', ...
        'doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', struct(...
                    'stimulusChromaticity', 'Mcone isolating',...
                    'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition));



    exampleLconeCenterRGCposition = [];
    exampleMconeCenterRGCposition = [];
    while (numel(exampleLconeCenterRGCposition) ~= 2)
        exampleLconeCenterRGCposition = input('Enter example L-cone center RGC position (e.g. [0.7 0.7]) :');
    end

    while (numel(exampleMconeCenterRGCposition) ~= 2)
        exampleMconeCenterRGCposition = input('Enter example M-cone center RGC position (e.g. [0.7 0.7]) :');
    end

    % Load the compute-ready mRGCMosaic
    load(fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

    % Analyze the surround cone mix for all mRGCs in theComputeReadyMRGCmosaic
    [surroundConeMix, theCenterMajorityConeType] = surroundConeMixForAllCellsInMosaic(...
            theComputeReadyMRGCmosaic, performSurroundAnalysisForConesExclusiveToTheSurround);
        
    % Determine visualized mRGCindices (only those whose surround cone mix is within the target range)
    mRGCindicesToVisualizeSTFsAcrossChromaticities = find(...
        (surroundConeMix >= targetRangeForSurroundConeMix(1)) & ...
        (surroundConeMix <= targetRangeForSurroundConeMix(2)));

    idxLconeCenter = find(theCenterMajorityConeType(mRGCindicesToVisualizeSTFsAcrossChromaticities) == cMosaic.LCONE_ID);
    idxMconeCenter = find(theCenterMajorityConeType(mRGCindicesToVisualizeSTFsAcrossChromaticities) == cMosaic.MCONE_ID);


    % Visualize the surround mix histograms for all L-center and
    % M-center cells in this mRGCmosaic
    pdfFileName = sprintf('SurroundConeMix_eccDegs_%2.1f_%2.1f.pdf', mosaicParams.eccDegs(1), mosaicParams.eccDegs(1));

     
    hFig = figure(555); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff, 'figPosition', [576 707]);
    plotTitle = '';
    MSreadyPlot.renderSurroundMixHistograms(theAxes{1,1}, ...
            surroundConeMix(find(theCenterMajorityConeType == cMosaic.LCONE_ID)), ...
            surroundConeMix(find(theCenterMajorityConeType == cMosaic.MCONE_ID)), ...
            plotTitle, ff, ...
            'targetRangeForSurroundConeMix', targetRangeForSurroundConeMix);


    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFig, 300);

    exportScaledFigureVersionForManuscript = true;
    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end


    % Visualize the locations of cells with a surround cone mix in the targetRangeForSurroundConeMix
    pdfFileNameLcenter = sprintf('SurroundConeMixOutlinedLcenterMRGClocationsWithinTargetRange_eccDegs_%2.1f_%2.1f.pdf', mosaicParams.eccDegs(1), mosaicParams.eccDegs(1));
    pdfFileNameMcenter = sprintf('SurroundConeMixOutlinedMcenterMRGClocationsWithinTargetRange_eccDegs_%2.1f_%2.1f.pdf', mosaicParams.eccDegs(1), mosaicParams.eccDegs(1));


    hFigLcenter = figure(556); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theLcenterAxes = MSreadyPlot.generateAxes(hFigLcenter,ff, 'figPosition', [268 127]);

    hFigMcenter = figure(557); clf;
    theMcenterAxes = MSreadyPlot.generateAxes(hFigMcenter,ff, 'figPosition', [969 127]);

    plotTitleLcenter = sprintf('L-center mRGCs with surround cone mix in [%2.2f-%2.2f]',targetRangeForSurroundConeMix(1), targetRangeForSurroundConeMix(2));
    plotTitleMcenter = sprintf('M-center mRGCs with surround cone mix in [%2.2f-%2.2f]',targetRangeForSurroundConeMix(1), targetRangeForSurroundConeMix(2));
    MSreadyPlot.renderIdentifiedMRGClocations(hFigLcenter, hFigMcenter, theLcenterAxes{1,1}, theMcenterAxes{1,1}, ...
            theComputeReadyMRGCmosaic, ...
            mRGCindicesToVisualizeSTFsAcrossChromaticities(idxLconeCenter), ...
            mRGCindicesToVisualizeSTFsAcrossChromaticities(idxMconeCenter), ...
            plotTitleLcenter, plotTitleMcenter, ff);


   
    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileNameLcenter);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigLcenter, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileNameLcenter);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileNameMcenter);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigMcenter, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileNameMcenter);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end


    MosaicPoolingOptimizer.contrastVisualSTFsAcrossDifferentChromaticities(...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName), ...
        fullfile(resourcesDirectory, mRGCMosaicAchromaticSTFresponsesFileName), ...
        fullfile(resourcesDirectory, mRGCMosaicLconeIsolatingSTFresponsesFileName), ...
        fullfile(resourcesDirectory, mRGCMosaicMconeIsolatingSTFresponsesFileName), ...
        opticsParamsForMRGCSTFs, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript,  ...
        'examinedRGCindices', mRGCindicesToVisualizeSTFsAcrossChromaticities, ...
        'tickSeparationArcMin', 2.0, ...
        'performSurroundAnalysisForConesExclusiveToTheSurround', performSurroundAnalysisForConesExclusiveToTheSurround, ...
        'generateVideoWithAllExaminedRGCs', false)

    
end


function [surroundConeMixForAllCells, theCenterMajorityConeTypeForAllCells] = surroundConeMixForAllCellsInMosaic(theMRGCmosaic, performSurroundAnalysisForConesExclusiveToTheSurround)

    surroundConeMixForAllCells = zeros(1, theMRGCmosaic.rgcsNum);
    theCenterMajorityConeTypeForAllCells = zeros(1, theMRGCmosaic.rgcsNum);

    parfor theRGCindex = 1:theMRGCmosaic.rgcsNum

        [theCenterMajorityConeType, ~,~, netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix] = MosaicPoolingOptimizer.analyzeCenterSurroundConeMix(...
         theMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround)

        surroundConeMixForAllCells(theRGCindex) = surroundConeMix;
        theCenterMajorityConeTypeForAllCells(theRGCindex) = theCenterMajorityConeType;
    end

end

function [visualizedConeIndices, theVisualizedConeXcoords, ...
          visualizedMRGCindices, theVisualizedMRGCXcoords, theROI] = extractVisualizedConeAndRGCindices(theMRGCMosaic, targetYdegs)

    % Define an ROI
    roiHeightDegs = 0.1;
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [theMRGCMosaic.eccentricityDegs(1) targetYdegs], ...
            'width', theMRGCMosaic.sizeDegs(1), ...
            'height', roiHeightDegs , ...
            'rotation', 0.0...
        ));

    visualizedConeIndices = theROI.indicesOfPointsInside(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs);
    visualizedMRGCindices = theROI.indicesOfPointsInside(theMRGCMosaic.rgcRFpositionsDegs);

    theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
    theVisualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));

    [~,idx] = sort(theVisualizedConeXcoords, 'ascend');
    visualizedConeIndices = visualizedConeIndices(idx);
    % Exclude S-cones
    idx = find(theMRGCMosaic.inputConeMosaic.coneTypes(visualizedConeIndices) == cMosaic.SCONE_ID);
    [~, idx] = setdiff(visualizedConeIndices, visualizedConeIndices(idx));
    visualizedConeIndices = visualizedConeIndices(idx);
    theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));

    [~,idx] = sort(theVisualizedMRGCXcoords, 'ascend');
    theVisualizedMRGCXcoords = theVisualizedMRGCXcoords(idx);
    visualizedMRGCindices = visualizedMRGCindices(idx);
end
