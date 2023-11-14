function contrastVisualSTFsAcrossDifferentChromaticities(...
        exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
        computeReadyMosaicFilename, ...
        mRGCMosaicAchromaticSTFresponsesFileName, ...
        mRGCMosaicLconeIsolatingSTFresponsesFileName, ...
        mRGCMosaicMconeIsolatingSTFresponsesFileName, ...
        opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('pdfDirectory', '', @ischar);
    p.addParameter('examinedRGCindices', [], @(x)(isnumeric(x)||(ischar(x))));
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.addParameter('performSurroundAnalysisForConesExclusiveToTheSurround', true, @islogical);
    p.addParameter('generateVideoWithAllExaminedRGCs', false, @islogical);

    p.parse(varargin{:});
    pdfDirectory = p.Results.pdfDirectory;
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

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

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


    plotSelectedRGCData(theComputeReadyMRGCmosaic, theLconeCenterData, theMconeCenterData, exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, opticsString)

end

function plotSelectedRGCData(theComputeReadyMRGCmosaic, theLconeCenterData, theMconeCenterData, exampleLconeCenterRGCposition, exampleMconeCenterRGCposition,  opticsString)

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
    exampleLconeCenterData = theLconeCenterData{idx};
    
    distances = sum((bsxfun(@minus, theMconeCenterRFpositions, exampleMconeCenterRGCposition)).^2,2);
    [~, idx] = min(distances(:));
    exampleMconeCenterData = theMconeCenterData{idx};


    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';

    % 1. Render the current population BPI plot for Lcenter RGCs
    hFigBPIpopulation = figure(10); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigBPIpopulation,ff);
    position = get(hFigBPIpopulation, 'OuterPosition');
    position(2) = 600;
    set(hFigBPIpopulation, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    MSreadyPlot.renderSTFBPIplot(ax, ...
        theLconeCenterBPIscatterData, ...
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
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,strrep(BPIpolulationPdfFileName, '.pdf', 'LcenterRGCs.pdf')), hFigBPIpopulation, 300);
    close(hFigBPIpopulation);

    % 1B. Render the current population BPI plot for Mcenter RGCs
    hFigBPIpopulation2 = figure(11); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigBPIpopulation2,ff);
    position = get(hFigBPIpopulation2, 'OuterPosition');
    position(2) = 600;
    set(hFigBPIpopulation2, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    MSreadyPlot.renderSTFBPIplot(ax, ...
        [], ...
        theMconeCenterBPIscatterData, ...
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
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,strrep(BPIpolulationPdfFileName, '.pdf', 'McenterRGCs.pdf')), hFigBPIpopulation2, 300);
    close(hFigBPIpopulation2);


    % 2A. Example L-cone RGC: STF magnitude spectra for the achromatic, L-cone
    % isolating and M-cone isolating gratings

    hFigMultiSTF = figure(20); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigMultiSTF,ff);
    set(hFigMultiSTF, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeMultiSTF(ax,  exampleLconeCenterData, ff);
    theSTFpdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'exampleLconeCenterRGC_STFs.pdf');
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot, theSTFpdfFileName), hFigMultiSTF, 300);
    close(hFigMultiSTF);

    % 2B. Example M-cone RGC: STF magnitude spectra for the achromatic, L-cone
    % isolating and M-cone isolating gratings

    hFigMultiSTF2 = figure(21); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigMultiSTF2,ff);
    set(hFigMultiSTF2, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeMultiSTF(ax,  exampleMconeCenterData, ff);
    theSTFpdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'exampleMconeCenterRGC_STFs.pdf');
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot, theSTFpdfFileName), hFigMultiSTF2, 300);
    close(hFigMultiSTF2);


    % 3A. Example L-cone RGC: Cone pooling map
    hFigConePoolingMap = figure(30); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigConePoolingMap,ff);
    set(hFigConePoolingMap, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, exampleLconeCenterData, cMosaic.LCONE_ID, ff);
    theConePoolingMapPdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'exampleLconeCenterRGC_ConePoolingMap.pdf');
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,theConePoolingMapPdfFileName), hFigConePoolingMap, 300);
    close(hFigConePoolingMap);

    % 3B. Example M-cone RGC: Cone pooling map
    hFigConePoolingMap2 = figure(31); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFigConePoolingMap2,ff);
    set(hFigConePoolingMap2, 'Color', [1 1 1], 'OuterPosition', position);
    ax = theAxes{1,1};

    visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, exampleMconeCenterData, cMosaic.MCONE_ID, ff);
    theConePoolingMapPdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'exampleMconeCenterRGC_ConePoolingMap.pdf');
    NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,theConePoolingMapPdfFileName), hFigConePoolingMap2, 300);
    close(hFigConePoolingMap2);

    % Generate paper-ready figures (scaled versions of the figures in nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end


function visualizeRetinalConePooling(ax, theComputeReadyMRGCmosaic, d, theCenterMajorityConeType, ff)
    theRearrangedAxes{1,1} = ax; % Plot the 2D cone pooling map
    theRearrangedAxes{1,2} = []; % Do not plot the X-line weighting function
    theRearrangedAxes{1,3} = []; % Do not plot the Y-line weighting function
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
     theSTFmagnitudeSpectra = theSTFmagnitudeSpectra / max(theSTFmagnitudeSpectra(:));
     theShadedSTFmagnitudeSpectrum = theSTFmagnitudeSpectra(:,1);
     theShadedSTFmagnitudeSpectrum = [];
     MSreadyPlot.renderMultiSTF(ax, d.spatialFrequenciesTested, ...
                 theSTFmagnitudeSpectra, theShadedSTFmagnitudeSpectrum, ...
                 theSTFphaseSpectra, theSTFColors, theLegends, sprintf('RGC #%d', d.theRGCindex), ff, ...
                 'noYLabel', false, ...
                 'noYTickLabel', true, ...
                 'visualizedSpatialFrequencyRange', [0.2 80]);
       
end
