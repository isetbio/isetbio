function contrastVisualSTFsAcrossDifferentChromaticities(...
        computeReadyMosaicFilename, ...
        mRGCMosaicAchromaticSTFresponsesFileName, ...
        mRGCMosaicLconeIsolatingSTFresponsesFileName, ...
        mRGCMosaicMconeIsolatingSTFresponsesFileName, ...
        opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('pdfDirectory', '', @ischar);
    p.addParameter('examinedRGCindices', [], @(x)(isnumeric(x)||(ischar(x))));
    p.addParameter('maxRGCsToIncludeWithinTheTargetRange', inf, @isscalar);
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.addParameter('performSurroundAnalysisForConesExclusiveToTheSurround', true, @islogical);
    p.addParameter('generateVideoWithAllExaminedRGCs', false, @islogical);

    p.parse(varargin{:});
    pdfDirectory = p.Results.pdfDirectory;
    examinedRGCindices = p.Results.examinedRGCindices;
    maxRGCsToIncludeWithinTheTargetRange = p.Results.maxRGCsToIncludeWithinTheTargetRange;
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

    %if (generateVideoWithAllExaminedRGCs)
    %    videoFileName = fullfile(pdfDirectory,sprintf('AchromaticLMconeIsolatingSTFanalysis_%s.mp4', opticsString));
    %    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    %    videoOBJ.FrameRate = 10;
    %    videoOBJ.Quality = 100;
    %    videoOBJ.open();
    %end

    analyzedLconeCenterRGCs = 0;
    analyzedMconeCenterRGCs = 0;

    for iRGC = 1:numel(examinedRGCindices)
        if (analyzedLconeCenterRGCs + analyzedMconeCenterRGCs > maxRGCsToIncludeWithinTheTargetRange)
            continue;
        end

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

    theLconeCenterBPIscatterData = [];
    theMconeCenterBPIscatterData = [];
    currentLconeCenterRGCindex = 0;
    currentMconeCenterRGCindex = 0;

    timeSupport = 0:(size(theMRGCMosaicAchromaticSTFresponses,3)-1);
    timeSupport = timeSupport/(numel(timeSupport)-1);
    timeSupportHR = timeSupport(1):0.01:timeSupport(end);
    temporalFrequency = 1.0;


    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';

    haveNotEncounteredLastLConeRGC = true;
    haveNotEncounteredLastMConeRGC = true;

    for iRGC = 1:numel(examinedRGCindices)

        % Only visualize the desired mRGCs
        %if (size(theLconeCenterBPIscatterData,1) + size(theMconeCenterBPIscatterData,1) == maxRGCsToIncludeWithinTheTargetRange) 

        % Only visualize the last L-cone center and the last M-cone center mRGCs
        onLastLconeCenterRGC = ((analyzedLconeCenterRGCs > 0) && (currentLconeCenterRGCindex == analyzedLconeCenterRGCs)) && (haveNotEncounteredLastLConeRGC);                             
        onLastMconeCenterRGC = ((analyzedMconeCenterRGCs > 0) && (currentMconeCenterRGCindex == analyzedMconeCenterRGCs)) && (haveNotEncounteredLastMConeRGC);

        if (onLastLconeCenterRGC || onLastMconeCenterRGC)
       % if (((analyzedLconeCenterRGCs > 0) && (currentLconeCenterRGCindex == analyzedLconeCenterRGCs)) && ((analyzedMconeCenterRGCs == 0)||(currentMconeCenterRGCindex ~= analyzedMconeCenterRGCs))) || ...
       %    (((analyzedMconeCenterRGCs > 0) && (currentMconeCenterRGCindex == analyzedMconeCenterRGCs)) && ((analyzedLconeCenterRGCs == 0)||(currentLconeCenterRGCindex ~= analyzedLconeCenterRGCs)))

            % 1. Render the current population BPI plot, outlining the BPI of
            % the last L-center and last M-center cell plotted
            hFigBPIpopulation = figure(10); clf;
            ff = MSreadyPlot.figureFormat('1x1 small');
            theAxes = MSreadyPlot.generateAxes(hFigBPIpopulation,ff);
            position = get(hFigBPIpopulation, 'OuterPosition');
            position(2) = 600;
            set(hFigBPIpopulation, 'Color', [1 1 1], 'OuterPosition', position);
            ax = theAxes{1,1};

            if (onLastLconeCenterRGC)
                haveNotEncounteredLastLConeRGC = false;
                MSreadyPlot.renderSTFBPIplot(ax, ...
                    theLconeCenterBPIscatterData, ...
                    [], ...
                    '', ff, ...
                    'gridless', false, ...
                    'opticsString', strrep(opticsString,'_', ':'), ...
                    'superimposeLeeShapleyData', true, ...
                    'showLegends', false, ...
                    'outlineLastPoint', true);
            else
                haveNotEncounteredLastMConeRGC = false;
                MSreadyPlot.renderSTFBPIplot(ax, ...
                    [], ...
                    theMconeCenterBPIscatterData, ...
                    '', ff, ...
                    'gridless', false, ...
                    'opticsString', strrep(opticsString,'_', ':'), ...
                    'superimposeLeeShapleyData', true, ...
                    'showLegends', false, ...
                    'outlineLastPoint', true);
            end

            
            % The pdf filename
            BPIpolulationPdfFileName = sprintf('BPIanalysisPopulation_eccDegs_%2.1f_%2.1f_%s.pdf', ...
               theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(2), strrep(strrep(strrep(opticsString, ', ', '_'), ':', '_'), ' ', ''));

            
            % 2. Render the STF magnitude spectra for the achromatic, L-cone
            % isolating and M-cone isolating gratings
            hFigMultiSTF = figure(20); clf;
            ff = MSreadyPlot.figureFormat('1x1 small');
            theAxes = MSreadyPlot.generateAxes(hFigMultiSTF,ff);
            set(hFigMultiSTF, 'Color', [1 1 1], 'OuterPosition', position);
            ax = theAxes{1,1};

            theSTFColors = [ ...
                0.2 0.2 0.2; ...
                1.0 0.0 0.0; ...
                0.0 0.8 0.0];

            % The STF magnitude spectra for the achromatic, L-cone isolating and
            % M-cone isolating gratings
            theSTFmagnitudeSpectra = [...
                theOptimalAchromaticSTFmagnitudeSpectrum(:) ...
                theOptimalLconeIsolatingSTFmagnitudeSpectrum(:) ...
                theOptimalMconeIsolatingSTFmagnitudeSpectrum(:)];

            % The STF phase spectra for the achromatic, L-cone isolating and
            % M-cone isolating gratings
            theSTFphaseSpectra = [...
                theOptimalAchromaticSTFphaseSpectrum(:) ...
                theOptimalLconeIsolatingSTFphaseSpectrum(:) ...
                theOptimalMconeIsolatingSTFphaseSpectrum(:)];

            theLegends = {'L+M+S', 'L-cone', 'M-cone'};
            theSTFmagnitudeSpectra = theSTFmagnitudeSpectra / max(theSTFmagnitudeSpectra(:));
            theShadedSTFmagnitudeSpectrum = theSTFmagnitudeSpectra(:,1);
            theShadedSTFmagnitudeSpectrum = [];
            MSreadyPlot.renderMultiSTF(ax, spatialFrequenciesTested, ...
                 theSTFmagnitudeSpectra, theShadedSTFmagnitudeSpectrum, ...
                 theSTFphaseSpectra, theSTFColors, theLegends, '', ff, ...
                 'noYLabel', false, ...
                 'noYTickLabel', true, ...
                 'visualizedSpatialFrequencyRange', [0.2 80]);
       
            lastfSTFpdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'LastVisualizedRGC_STFs.pdf');


            % 3. Render the cone pooling map
            hFigConePoolingMap = figure(30); clf;
            ff = MSreadyPlot.figureFormat('1x1 small');
            theAxes = MSreadyPlot.generateAxes(hFigConePoolingMap,ff);
            set(hFigConePoolingMap, 'Color', [1 1 1], 'OuterPosition', position);
            ax = theAxes{1,1};

            theRearrangedAxes{1,1} = ax; % Plot the 2D cone pooling map
            theRearrangedAxes{1,2} = []; % Do not plot the X-line weighting function
            theRearrangedAxes{1,3} = []; % Do not plot the Y-line weighting function
            tickSeparationArcMin = 3;

            theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapOfRGCwithIndex(...
                theRGCindex, ...
                'theAxes', theRearrangedAxes, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'normalizedPeakSurroundSensitivity', 0.4, ...
                'gridlessLineWeightingFunctions', true, ...
                'withFigureFormat', ff, ...
                'regenerateVisualizationCache', false);

            centerL = netCenterLconeWeight / (netCenterLconeWeight+netCenterMconeWeight);
            centerM = netCenterMconeWeight / (netCenterLconeWeight+netCenterMconeWeight);
            surroundL = netSurroundLconeWeight / (netSurroundLconeWeight+netSurroundMconeWeight);
            surroundM = netSurroundMconeWeight / (netSurroundLconeWeight+netSurroundMconeWeight);
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

            lastConePoolingMapPdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'LastVisualizedRGC_ConePoolingMap.pdf');


            % 4. Render the X-line weighting function
            hFigLineWeightingFunctionX = figure(40); clf;
            ff = MSreadyPlot.figureFormat('1x1 small');
            theAxes = MSreadyPlot.generateAxes(hFigLineWeightingFunctionX,ff);
            set(hFigLineWeightingFunctionX, 'Color', [1 1 1], 'OuterPosition', position);
            ax = theAxes{1,1};

            theRearrangedAxes{1,1} = ax;
            theRearrangedAxes{1,2} = ax;
            theRearrangedAxes{1,3} = []; % Do not plot the Y-line weighting function

            theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapOfRGCwithIndex(...
                theRGCindex, ...
                'theAxes', theRearrangedAxes, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'normalizedPeakSurroundSensitivity', 0.4, ...
                'gridlessLineWeightingFunctions', false, ...
                'noYLabelLineWeightingFunctions', false, ...
                'withFigureFormat', ff, ...
                'regenerateVisualizationCache', false);

            title(theRearrangedAxes{1,1}, sprintf('RGC #%d', theRGCindex),'FontSize', ff.titleFontSize, ...
                'FontWeight', ff.titleFontWeight, 'Color', titleColor);

            lastConePoolingLineWeightingFunctionXPdfFileName = strrep(BPIpolulationPdfFileName, '.pdf', 'LastVisualizedRGC_ConeLineWeightingFunctionX.pdf');

            % Print
            NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,BPIpolulationPdfFileName), hFigBPIpopulation, 300);
            NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,lastfSTFpdfFileName), hFigMultiSTF, 300);
            NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,lastConePoolingMapPdfFileName), hFigConePoolingMap, 300);
            NicePlot.exportFigToPDF(fullfile(rawFiguresRoot,lastConePoolingLineWeightingFunctionXPdfFileName), hFigLineWeightingFunctionX, 300);

            % Generate paper-ready figures (scaled versions of the figures i
            % nrawFiguresRoot directory) which are stored in the PaperReady folder
            PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
            commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
            system(commandString);

            disp('Hit enter to continue');
            pause
        end
       
        allDoneWithLconeCenterRGCs = (analyzedLconeCenterRGCs == 0) || (analyzedLconeCenterRGCs > 0) && (haveNotEncounteredLastLConeRGC == false);
        allDoneWithMconeCenterRGCs = (analyzedMconeCenterRGCs == 0) || (analyzedMconeCenterRGCs > 0) && (haveNotEncounteredLastMConeRGC == false);
        if (allDoneWithLconeCenterRGCs && allDoneWithMconeCenterRGCs)
            continue;
        end



        % Get theRGCindex
        theRGCindex = examinedRGCindices(iRGC);

        % Get the achromatic STFs across all orientations for this RGC
        [~, ~, theHighestExtensionOrientation, ...
            theOptimalAchromaticSTFmagnitudeSpectrum, theOptimalAchromaticSTFphaseSpectrum] = MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
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

        for iSF = 1:numel(spatialFrequenciesTested)
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
         netSurroundLconeWeight, netSurroundMconeWeight] = analyzeCenterSurroundConeMix(...
            theComputeReadyMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround);

        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            currentLconeCenterRGCindex = currentLconeCenterRGCindex + 1;
            theLconeCenterBPIscatterData(currentLconeCenterRGCindex,:) = [theAchromaticSTFBandPassIndex theLconeIsolatingSTFBandPassIndex];
        else
            currentMconeCenterRGCindex = currentMconeCenterRGCindex + 1;
            theMconeCenterBPIscatterData(currentMconeCenterRGCindex+1,:) = [theAchromaticSTFBandPassIndex theMconeIsolatingSTFBandPassIndex];
        end
    
        %if (generateVideoWithAllExaminedRGCs)
        %    videoOBJ.writeVideo(getframe(hFig));
        %end

    end % iRGC

    %if (generateVideoWithAllExaminedRGCs)
    %    videoOBJ.close();
    %end

    

end

function [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight] = analyzeCenterSurroundConeMix(...
         theMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround)

    [theCenterConeTypeNetWeights, ~, theCenterMajorityConeType, ...
        theCenterConeTypes, theCenterConeIndices] = theMRGCmosaic.centerConeTypeWeights(theRGCindex);
    netCenterLconeWeight = theCenterConeTypeNetWeights(find(theCenterConeTypes == cMosaic.LCONE_ID));
    netCenterMconeWeight = theCenterConeTypeNetWeights(find(theCenterConeTypes == cMosaic.MCONE_ID));

    if (performSurroundAnalysisForConesExclusiveToTheSurround)
        [~, theSurroundConeTypeNetWeights] = theMRGCmosaic.surroundConeTypeWeights(theRGCindex, theCenterConeIndices);
    else
        theSurroundConeTypeNetWeights = theMRGCmosaic.surroundConeTypeWeights(theRGCindex, []);
    end

    netSurroundLconeWeight = theSurroundConeTypeNetWeights(cMosaic.LCONE_ID);
    netSurroundMconeWeight = theSurroundConeTypeNetWeights(cMosaic.MCONE_ID);
end

