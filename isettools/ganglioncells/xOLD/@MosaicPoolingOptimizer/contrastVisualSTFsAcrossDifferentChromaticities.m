function contrastVisualSTFsAcrossDifferentChromaticities(...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        theComputeReadyMRGCmosaic, ...
        mRGCMosaicAchromaticSTFresponsesFileName, ...
        mRGCMosaicLconeIsolatingSTFresponsesFileName, ...
        mRGCMosaicMconeIsolatingSTFresponsesFileName, ...
        opticsParams, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript, ...
        varargin)

    % Parse input
    p = inputParser;
    p.addParameter('examinedRGCindices', [], @(x)(isnumeric(x)||(ischar(x))));
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.addParameter('performSurroundAnalysisForConesExclusiveToTheSurround', true, @islogical);
    p.addParameter('generateVideoWithAllExaminedRGCs', false, @islogical);

    p.parse(varargin{:});
    examinedRGCindices = p.Results.examinedRGCindices;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    performSurroundAnalysisForConesExclusiveToTheSurround = p.Results.performSurroundAnalysisForConesExclusiveToTheSurround;
    generateVideoWithAllExaminedRGCs = p.Results.generateVideoWithAllExaminedRGCs;

    % Load the achromatic STF responses
    load(mRGCMosaicAchromaticSTFresponsesFileName, 'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
        'orientationsTested', 'spatialFrequenciesTested');
    theMRGCMosaicAchromaticSTFresponses = theMRGCMosaicSTFresponses;
    clear 'theMRGCMosaicSTFresponses';

    % Load the L-cone isolating STF responses
    load(mRGCMosaicLconeIsolatingSTFresponsesFileName, 'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
        'orientationsTested', 'spatialFrequenciesTested');
    theMRGCMosaicLconeIsolatingSTFresponses = theMRGCMosaicSTFresponses;
    clear 'theMRGCMosaicSTFresponses';


    % Load the M-cone isolating STF responses
    load(mRGCMosaicMconeIsolatingSTFresponsesFileName, 'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
        'orientationsTested', 'spatialFrequenciesTested');
    theMRGCMosaicMconeIsolatingSTFresponses= theMRGCMosaicSTFresponses;
    clear 'theMRGCMosaicSTFresponses';

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

    
    % Generate the mRGCMosaic visualization cache
    xSupport = [];
    ySupport = []; 
    centerSubregionContourSamples = 32;
    contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';
    theComputeReadyMRGCmosaic.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);
  

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

    timeSupport = 0:(size(theMRGCMosaicAchromaticSTFresponses,3)-1);
    timeSupport = timeSupport/(numel(timeSupport)-1);
    timeSupportHR = timeSupport(1):0.01:timeSupport(end);
    temporalFrequency = 1.0;

    theLconeCenterData = cell(1, analyzedLconeCenterRGCs);
    theMconeCenterData = cell(1, analyzedMconeCenterRGCs);

    currentLconeCenterRGCindex = 0;
    currentMconeCenterRGCindex = 0;

    for iRGC = 1:numel(examinedRGCindices)

        fprintf('Computing L-, M-cone isolating and achromatic STFs for cell %d of %d\n', iRGC, numel(examinedRGCindices));

        % Get theRGCindex
        theRGCindex = examinedRGCindices(iRGC);
        
        % Get the achromatic STFs across all orientations for this RGC
        [~, ~, theHighestExtensionOrientation, theOptimalAchromaticSTFmagnitudeSpectrum, theOptimalAchromaticSTFphaseSpectrum] = ...
            MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
                    orientationsTested, spatialFrequenciesTested, ...
                    squeeze(theMRGCMosaicAchromaticSTFresponses(:,:,:, theRGCindex)));
        % Determine the orientation index that resulted in the optimal achromatic STF
        iOri = find(orientationsTested == theHighestExtensionOrientation);

        % Compute the L- and M-cone isolating STFs for the orientation of
        % the optimal achromatic STF
        theOptimalLconeIsolatingSTFmagnitudeSpectrum = zeros(1, numel(spatialFrequenciesTested));
        theOptimalMconeIsolatingSTFmagnitudeSpectrum = zeros(1, numel(spatialFrequenciesTested));
        theOptimalLconeIsolatingSTFphaseSpectrum = zeros(1, numel(spatialFrequenciesTested));
        theOptimalMconeIsolatingSTFphaseSpectrum  = zeros(1, numel(spatialFrequenciesTested));

        parfor iSF = 1:numel(spatialFrequenciesTested)
            % Retrieve L-cone isolating responses for all frames 
            theResponseTimeSeries = squeeze(theMRGCMosaicLconeIsolatingSTFresponses(iOri, iSF, :, theRGCindex));
            % Fit a sinusoid to extract the magnitude and phase
            [~, fittedParams] = MosaicPoolingOptimizer.fitSinusoidToResponseTimeSeries(...
                timeSupport, theResponseTimeSeries, temporalFrequency, timeSupportHR);
            theOptimalLconeIsolatingSTFmagnitudeSpectrum(iSF) = fittedParams(1);
            theOptimalLconeIsolatingSTFphaseSpectrum(iSF) = fittedParams(2);

            % Retrieve M-cone isolating responses for all frames 
            theResponseTimeSeries = squeeze(theMRGCMosaicMconeIsolatingSTFresponses(iOri, iSF, :, theRGCindex));
            % Fit a sinusoid to extract the magnitude and phase
            [~, fittedParams] = MosaicPoolingOptimizer.fitSinusoidToResponseTimeSeries(...
                timeSupport, theResponseTimeSeries, temporalFrequency, timeSupportHR);
            theOptimalMconeIsolatingSTFmagnitudeSpectrum(iSF) = fittedParams(1);
            theOptimalMconeIsolatingSTFphaseSpectrum(iSF) = fittedParams(2);
        end

        % Apply correction factors to account for low L/M cone isolating stimulus contrast
        % compared to the achromatic stimulus contrast
        theOptimalLconeIsolatingSTFmagnitudeSpectrum = theOptimalLconeIsolatingSTFmagnitudeSpectrum * LconeStimulusContrastCorrection;
        theOptimalMconeIsolatingSTFmagnitudeSpectrum = theOptimalMconeIsolatingSTFmagnitudeSpectrum * MconeStimulusContrastCorrection;

        theLconeIsolatingSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theOptimalLconeIsolatingSTFmagnitudeSpectrum);
        theMconeIsolatingSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theOptimalMconeIsolatingSTFmagnitudeSpectrum);
        theAchromaticSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theOptimalAchromaticSTFmagnitudeSpectrum);

        % Analyze cone weights to center and surround
        [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix] = MosaicPoolingOptimizer.analyzeCenterSurroundConeMix(...
            theComputeReadyMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround);

        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            currentLconeCenterRGCindex = currentLconeCenterRGCindex + 1;
            theLconeCenterData{currentLconeCenterRGCindex} = struct(...
                'theRGCindex', theRGCindex, ...
                'theHighestExtensionOrientation', theHighestExtensionOrientation, ...
                'bandpassIndices', [theAchromaticSTFBandPassIndex theLconeIsolatingSTFBandPassIndex], ...
                'netCenterLconeWeight', netCenterLconeWeight, ...
                'netCenterMconeWeight', netCenterMconeWeight, ...
                'netSurroundLconeWeight', netSurroundLconeWeight, ...
                'netSurroundMconeWeight', netSurroundMconeWeight, ...
                'surroundConeMix',  surroundConeMix, ...
                'spatialFrequenciesTested', spatialFrequenciesTested, ...
                'achromaticSTFamplitude', theOptimalAchromaticSTFmagnitudeSpectrum, ...
                'LconeIsolatingSTFamplitude', theOptimalLconeIsolatingSTFmagnitudeSpectrum, ...
                'MconeIsolatingSTFmamplitude', theOptimalMconeIsolatingSTFmagnitudeSpectrum, ...
                'achromaticSTFphase', theOptimalAchromaticSTFphaseSpectrum, ...
                'LconeIsolatingSTFphase', theOptimalLconeIsolatingSTFphaseSpectrum, ...
                'MconeIsolatingSTFphase', theOptimalMconeIsolatingSTFphaseSpectrum);

        else
            currentMconeCenterRGCindex = currentMconeCenterRGCindex + 1;
            theMconeCenterData{currentMconeCenterRGCindex} = struct(...
                'theRGCindex', theRGCindex, ...
                'theHighestExtensionOrientation', theHighestExtensionOrientation, ...
                'bandpassIndices', [theAchromaticSTFBandPassIndex theMconeIsolatingSTFBandPassIndex], ...
                'netCenterLconeWeight', netCenterLconeWeight, ...
                'netCenterMconeWeight', netCenterMconeWeight, ...
                'netSurroundLconeWeight', netSurroundLconeWeight, ...
                'netSurroundMconeWeight', netSurroundMconeWeight, ...
                'surroundConeMix',  surroundConeMix, ...
                'spatialFrequenciesTested', spatialFrequenciesTested, ...
                'achromaticSTFamplitude', theOptimalAchromaticSTFmagnitudeSpectrum, ...
                'LconeIsolatingSTFamplitude', theOptimalLconeIsolatingSTFmagnitudeSpectrum, ...
                'MconeIsolatingSTFmamplitude', theOptimalMconeIsolatingSTFmagnitudeSpectrum, ...
                'achromaticSTFphase', theOptimalAchromaticSTFphaseSpectrum, ...
                'LconeIsolatingSTFphase', theOptimalLconeIsolatingSTFphaseSpectrum, ...
                'MconeIsolatingSTFphase', theOptimalMconeIsolatingSTFphaseSpectrum);
        end
    end % iRGC


    plotSelectedRGCData(theComputeReadyMRGCmosaic, ...
        theLconeCenterData, theMconeCenterData, ...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

end

function plotSelectedRGCData(theComputeReadyMRGCmosaic, theLconeCenterData, ...
    theMconeCenterData, exampleLconeCenterRGCposition, ...
    exampleMconeCenterRGCposition,  opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    lConeCenterRGCsNum = numel(theLconeCenterData);
    mConeCenterRGCsNum = numel(theMconeCenterData);
    theLconeCenterBPIscatterData = zeros(lConeCenterRGCsNum ,2);
    theMconeCenterBPIscatterData = zeros(mConeCenterRGCsNum ,2);

    theLconeCenterRFpositions = zeros(lConeCenterRGCsNum ,2);
    theMconeCenterRFpositions = zeros(mConeCenterRGCsNum ,2);

    for iLconeRGC = 1:lConeCenterRGCsNum
        d = theLconeCenterData{iLconeRGC};
        theLconeCenterBPIscatterData(iLconeRGC,:) = d.bandpassIndices;
        theLconeCenterRFpositions(iLconeRGC,:) = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(d.theRGCindex,:);
    end

    for iMconeRGC = 1:mConeCenterRGCsNum
        d = theMconeCenterData{iMconeRGC};
        theMconeCenterBPIscatterData(iMconeRGC,:) = d.bandpassIndices;
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


    % ============== Export Figures ==========================

    % 1. Render the current population BPI plot for Lcenter RGCs
    renderBPIfigure(10, theLconeCenterBPIscatterData, 'L', opticsString, theComputeReadyMRGCmosaic,  ...
        rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);
    renderBPIfigure(11, theMconeCenterBPIscatterData, 'M', opticsString, theComputeReadyMRGCmosaic,  ...
        rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    % 2. STF magnitude spectra for the achromatic, L-cone isolating and M-cone isolating gratings
    if (~isempty(exampleLconeCenterData))
        renderSTFsFigure(20, theComputeReadyMRGCmosaic, exampleLconeCenterData, 'L', ...
            opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);
    end
    if (~isempty(exampleMconeCenterData))
        renderSTFsFigure(21, theComputeReadyMRGCmosaic, exampleMconeCenterData, 'M', ...
            opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);
    end


    % 3A. Example L-cone RGC: Cone pooling 2D map
    if (~isempty(exampleLconeCenterData))
        renderConePoolingMapFigure(30, theComputeReadyMRGCmosaic, exampleLconeCenterData, cMosaic.LCONE_ID, 'L', ...
            opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);

    end
    if (~isempty(exampleMconeCenterData))
        renderConePoolingMapFigure(31, theComputeReadyMRGCmosaic, exampleMconeCenterData, cMosaic.MCONE_ID, 'M', ...
            opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);

    end

end


function renderBPIfigure(figNo, theBPIscatterData, centerConeType, ...
    opticsString, theComputeReadyMRGCmosaic, ...
    rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    hFigBPIpopulation = figure(figNo); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigBPIpopulation,ff);
    position = get(hFigBPIpopulation, 'OuterPosition');
    position(2) = 600;
    set(hFigBPIpopulation, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    MSreadyPlot.renderSTFBPIplot(ax, ...
        theBPIscatterData, ...
        [], ...
        '', ff, ...
        'gridless', false, ...
        'opticsString', strrep(opticsString,'_', ':'), ...
        'superimposeLeeShapleyData', true, ...
        'showLegends', false, ...
        'outlineLastPoint', true);
    
    % The pdf filename
    BPIpolulationPdfFileName = sprintf('BPIanalysisPopulation_eccDegs_%2.1f_%2.1f%s.pdf', ...
       theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(2), strrep(strrep(strrep(opticsString, ', ', '_'), ':', '_'), ' ', ''));

    % Print
    theConeBPIfigureName = strrep(BPIpolulationPdfFileName, '.pdf', sprintf('%scenterRGCs.pdf', centerConeType));
    pdfFileNameUnscaled = fullfile(rawFiguresRoot, theConeBPIfigureName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigBPIpopulation, 300);
    
    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, theConeBPIfigureName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end
    close(hFigBPIpopulation);
end



function renderSTFsFigure(figNo, theComputeReadyMRGCmosaic, exampleConeCenterData, coneType, ...
    opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    hFigMultiSTF = figure(figNo); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigMultiSTF,ff);
    position = get(hFigMultiSTF, 'OuterPosition');
    position(2) = 100;
    set(hFigMultiSTF, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeMultiSTF(ax,  exampleConeCenterData, ff);

    % The pdf filename
    BPIpolulationPdfFileName = sprintf('BPIanalysisPopulation_eccDegs_%2.1f_%2.1f%s.pdf', ...
       theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(2), strrep(strrep(strrep(opticsString, ', ', '_'), ':', '_'), ' ', ''));

    % Print
    pdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', sprintf('example%sconeCenterRGC_STFs.pdf', coneType));
    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled,  hFigMultiSTF, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    close(hFigMultiSTF);
end



function renderConePoolingMapFigure(figNo, theComputeReadyMRGCmosaic, ...
    exampleConeCenterData, coneID, coneType, ...
    opticsString, rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    % The 2D cone pooling map
    hFigConePoolingMap = figure(figNo); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigConePoolingMap,ff);
    
    position = get(hFigConePoolingMap, 'OuterPosition');
    position(2) = 100;

    set(hFigConePoolingMap, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, exampleConeCenterData, coneID, ff, false);

     % The pdf filename
    BPIpolulationPdfFileName = sprintf('BPIanalysisPopulation_eccDegs_%2.1f_%2.1f%s.pdf', ...
       theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(2), strrep(strrep(strrep(opticsString, ', ', '_'), ':', '_'), ' ', ''));

    % Print
    pdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', sprintf('example%sconeCenterRGC_ConePoolingMap.pdf', coneType));
    pdfFileNameUnscaled = fullfile(rawFiguresRoot,pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigConePoolingMap, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    close(hFigConePoolingMap);

    % The 1D line weighting function
    hFigConePoolingMap2 = figure(40); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigConePoolingMap2,ff);

    position = get(hFigConePoolingMap2, 'OuterPosition');
    position(2) = 400;

    set(hFigConePoolingMap2, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, exampleConeCenterData, coneID, ff, true);

    % Print
    pdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', sprintf('example%sconeCenterRGC_ConePoolingLineWeightingFunction.pdf', coneType));
    pdfFileNameUnscaled = fullfile(rawFiguresRoot,pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigConePoolingMap2, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    close(hFigConePoolingMap2);
end


function visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, d, theCenterMajorityConeType, ff, onlyVisualizeTheXprofile)

    theRearrangedAxes{1,1} = ax; % Plot the 2D cone pooling map
  
    if (onlyVisualizeTheXprofile)
        theRearrangedAxes{1,2} = ax; % Plot the X-line weighting function
        theRearrangedAxes{1,3} = []; % Do not plot the Y-line weighting function
    else
        theRearrangedAxes{1,2} = []; % Do not plot the X-line weighting function
        theRearrangedAxes{1,3} = []; % Do not plot the Y-line weighting function
    end

    tickSeparationArcMin = 3;

    theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapOfRGCwithIndex(...
        d.theRGCindex, ...
        'theAxes', theRearrangedAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'surroundSaturationLevel', 0.1, ...
        'normalizedPeakSurroundSensitivity', 0.4, ...
        'gridlessLineWeightingFunctions', true, ...
        'withFigureFormat', ff, ...
        'regenerateVisualizationCache', false);

     centerL = d.netCenterLconeWeight / (d.netCenterLconeWeight + d.netCenterMconeWeight);
     centerM = d.netCenterMconeWeight / (d.netCenterLconeWeight + d.netCenterMconeWeight);
     surroundL = d.netSurroundLconeWeight / (d.netSurroundLconeWeight + d.netSurroundMconeWeight);
     surroundM = d.netSurroundMconeWeight / (d.netSurroundLconeWeight + d.netSurroundMconeWeight);
     if ((centerL == 0) || (centerM == 0))
        if (centerL == 0)
            coneWeightsInfo = sprintf('L/M ratio: 0/1 (cntr), %0.2f/%0.2f (srrnd)', ...
                surroundL, surroundM);
        else
            coneWeightsInfo = sprintf('L/M ratio: 1/0 (cntr), %0.2f/%0.2f (srrnd)', ...
                surroundL, surroundM);
        end
     else
        coneWeightsInfo = sprintf('Lw/Mm: %0.1f/%0.1f (cntr), %0.2f/%0.2f (srrnd)', ...
            centerL, centerM, surroundL, surroundM);
     end
    
    if (theCenterMajorityConeType == cMosaic.LCONE_ID)
        titleColor = [0.8 0 0];
    else
        titleColor = [0 0.7 0.0];
    end
    title(theRearrangedAxes{1,1}, coneWeightsInfo,'FontSize', ff.titleFontSize, ...
                'FontWeight', ff.titleFontWeight, 'Color', titleColor);
end




function visualizeMultiSTF(ax,  d, ff)
    theSTFColors = [ ...
        0.2 0.2 0.2; ...
        1.0 0.0 0.0; ...
        0.0 0.8 0.0];

    % The STF magnitude spectra for the achromatic, L-cone isolating and
    % M-cone isolating gratings
    theSTFmagnitudeSpectra = [...
         d.achromaticSTFamplitude(:) ...
         d.LconeIsolatingSTFamplitude(:) ...
         d.MconeIsolatingSTFmamplitude(:)];

     % The STF phase spectra for the achromatic, L-cone isolating and
     % M-cone isolating gratings
     theSTFphaseSpectra = [...
         d.achromaticSTFphase(:) ...
         d.LconeIsolatingSTFphase(:) ...
         d.MconeIsolatingSTFphase(:)];

     theLegends = {'L+M+S', 'L-cone', 'M-cone'};
     theShadedSTFmagnitudeSpectrum = theSTFmagnitudeSpectra(:,1);
     theShadedSTFmagnitudeSpectrum = [];
     MSreadyPlot.renderMultiSTF(ax, d.spatialFrequenciesTested, ...
                 theSTFmagnitudeSpectra, theShadedSTFmagnitudeSpectrum, ...
                 theSTFphaseSpectra, theSTFColors, theLegends, sprintf('RGC #%d (ori: %d)', d.theRGCindex, d.theHighestExtensionOrientation), ff, ...
                 'noYLabel', false, ...
                 'noYTickLabel', false, ...
                 'visualizedSpatialFrequencyRange', [0.2 80]);
       
end
