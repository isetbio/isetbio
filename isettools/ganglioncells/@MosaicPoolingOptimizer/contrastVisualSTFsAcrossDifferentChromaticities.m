function contrastVisualSTFsAcrossDifferentChromaticities(...
        computeReadyMosaicFilename, ...
        mRGCMosaicAchromaticSTFresponsesFileName, ...
        mRGCMosaicLconeIsolatingSTFresponsesFileName, ...
        mRGCMosaicMconeIsolatingSTFresponsesFileName, ...
        opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('pdfDirectory', '', @ischar);
    p.addParameter('visualizedRGCindices', [], @(x)(isnumeric(x)||(ischar(x))));
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.parse(varargin{:});
    pdfDirectory = p.Results.pdfDirectory;
    visualizedRGCindices = p.Results.visualizedRGCindices;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;


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


    % achromaticStimulusRMSContrast = achromaticStimulusContrast * norm(achromaticStimulusConeContrasts(1:2));
    LconeStimulusContrastCorrection =  achromaticStimulusContrast * achromaticStimulusConeContrasts(1) / (LconeIsolatingStimulusContrast * LconeIsolatingStimulusConeContrasts(1));
    MconeStimulusContrastCorrection =  achromaticStimulusContrast * achromaticStimulusConeContrasts(2) / (MconeIsolatingStimulusContrast * MconeIsolatingStimulusConeContrasts(2));

    
    
    % Generate the visualization cache
    xSupport = [];
    ySupport = []; 
    centerSubregionContourSamples = 32;
    contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';
    theComputeReadyMRGCmosaic.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);
  

    if (isempty(visualizedRGCindices))
        visualizedRGCindices = 1:theComputeReadyMRGCmosaic.rgcsNum;
    end

    if (opticsParams.examinedSubjectRankOrder == 0)
       opticsString = sprintf('diffr-limited, pupil_%2.0fmm, defocus_%1.2fD', opticsParams.pupilDiameterMM, opticsParams.refractiveErrorDiopters);
    else
       opticsString = sprintf('physio-optics, defocus_%1.2fD', opticsParams.refractiveErrorDiopters);
    end

    videoFileName = fullfile(pdfDirectory,sprintf('AchromaticLMconeIsolatingSTFanalysis_%s.mp4', opticsString));
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    theLconeCenterBPIscatterData = [];
    theMconeCenterBPIscatterData = [];

    for iRGC = 1:numel(visualizedRGCindices)

        % Get theRGCindex
        theRGCindex = visualizedRGCindices(iRGC);

        % Get all the achromatic STFs for this RGC
        [~, ~, theHighestExtensionOrientation, theUnscaledOptimalAchromaticSTF] = MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
                    orientationsTested, spatialFrequenciesTested, ...
                    squeeze(theMRGCMosaicAchromaticSTFresponses(:,:,:, theRGCindex)));
        % Determine the orientation index that resulted in the optimal achromatic STF
        iOri = find(orientationsTested == theHighestExtensionOrientation);

        % Compute the optimal L- and M-cone isolating STFs
        theUnscaledOptimalLconeIsolatingSTF = zeros(1, numel(spatialFrequenciesTested));
        theUnscaledOptimalMconeIsolatingSTF = zeros(1, numel(spatialFrequenciesTested));

        for iSF = 1:numel(spatialFrequenciesTested)
            % Retrieve L-cone isolating responses for all frames 
            theResponseTimeSeries = squeeze(theMRGCMosaicLconeIsolatingSTFresponses(iOri, iSF, :, theRGCindex));
            % STF as the amplitude of modulation of the response 0.5*(max-min)
            theUnscaledOptimalLconeIsolatingSTF(iSF) = 0.5*(max(theResponseTimeSeries(:)) - min(theResponseTimeSeries(:)));
            % Retrieve M-cone isolating responses for all frames 
            theResponseTimeSeries = squeeze(theMRGCMosaicMconeIsolatingSTFresponses(iOri, iSF, :, theRGCindex));
            % STF as the amplitude of modulation of the response (max-min)
            theUnscaledOptimalMconeIsolatingSTF(iSF) = 0.5*(max(theResponseTimeSeries(:)) - min(theResponseTimeSeries(:)));
        end

        % Correct for stimulus contrast
        theUnscaledOptimalLconeIsolatingSTF = theUnscaledOptimalLconeIsolatingSTF * LconeStimulusContrastCorrection;
        theUnscaledOptimalMconeIsolatingSTF = theUnscaledOptimalMconeIsolatingSTF * MconeStimulusContrastCorrection;

        theLconeIsolatingSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theUnscaledOptimalLconeIsolatingSTF);
        theMconeIsolatingSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theUnscaledOptimalMconeIsolatingSTF);
        theAchromaticSTFBandPassIndex = MosaicPoolingOptimizer.STFbandpassIndex(spatialFrequenciesTested, theUnscaledOptimalAchromaticSTF);

        % Analyze cone weights to center and surround
        [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight] = ...
            analyzeSurroundConeMixing(theComputeReadyMRGCmosaic, theRGCindex);

        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            theLconeCenterBPIscatterData(size(theLconeCenterBPIscatterData,1)+1,:) = [theAchromaticSTFBandPassIndex theLconeIsolatingSTFBandPassIndex];
        else
            theMconeCenterBPIscatterData(size(theMconeCenterBPIscatterData,1)+1,:) = [theAchromaticSTFBandPassIndex theMconeIsolatingSTFBandPassIndex];
        end

        % Render the current BPI plot
        hFig = figure(2); clf;
        ff = MSreadyPlot.figureFormat('1x1 small');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        position = get(hFig, 'OuterPosition');
        position(2) = 600;
        set(hFig, 'Color', [1 1 1], 'OuterPosition', position);
        ax = theAxes{1,1};
        MSreadyPlot.renderSTFBPIplot(ax, ...
            theLconeCenterBPIscatterData, ...
            theMconeCenterBPIscatterData, ...
            '', ff, ...
            'gridless', false, ...
            'opticsString', strrep(opticsString,'_', ':'), ...
            'superimposeLeeShapleyData', true);

        BPIanalysisFileNamePDF = sprintf('BPIanalysis_%s.pdf', strrep(strrep(opticsString, ', ', '_'), ':', '_'));
        NicePlot.exportFigToPDF(fullfile(pdfDirectory,BPIanalysisFileNamePDF), hFig, 300);

        theSTFs = [theUnscaledOptimalAchromaticSTF(:) theUnscaledOptimalLconeIsolatingSTF(:) theUnscaledOptimalMconeIsolatingSTF(:)];
        theNormalizedSTFs = theSTFs / max(theSTFs(:));
        theLegends = {'L+M+S', 'L-cone', 'M-cone'};

        hFig = figure(1); clf;
        ff = MSreadyPlot.figureFormat('1x4 RF poster');
        ff.titleFontSize = 16;
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);
        axPSF = theAxes{1,3};
        axAchromaticLMisolatingSTFs = theAxes{1,4};

        theRearrangedAxes = theAxes;
        theRearrangedAxes{1,3} = [];

        

        theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapOfRGCwithIndex(...
            theRGCindex, ...
            'theAxes', theRearrangedAxes, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', 0.4, ...
            'gridlessLineWeightingFuncions', true, ...
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


        % Visualize the vLambda weighted PSF
        MosaicPoolingOptimizer.visualizeVlambdaWeightedPSF(...
            theComputeReadyMRGCmosaic, opticsParams, ...
            'axesHandle', axPSF, ...
            'figureFormat', ff);


        theSTFColors = [ ...
            0.2 0.2 0.2; ...
            1.0 0.0 0.0; ...
            0.0 0.8 0.0];


        MSreadyPlot.renderMultiSTF(axAchromaticLMisolatingSTFs, spatialFrequenciesTested, ...
                 theNormalizedSTFs, theNormalizedSTFs(:,1), theSTFColors, theLegends, strrep(opticsString, '_', ' '), ff, ...
                 'noYLabel', false, ...
                 'noYTickLabel', true, ...
                 'visualizedSpatialFrequencyRange', [min(spatialFrequenciesTested) max(spatialFrequenciesTested)]);
        

        if (theCenterMajorityConeType == cMosaic.LCONE_ID)
            fName = sprintf('AchromaticLMconeIsolatingSTFanalysis_RGC_%d_LconeCenter_%s.pdf', theRGCindex, strrep(strrep(opticsString, ', ', '_'), ':', '_'));
        elseif (theCenterMajorityConeType == cMosaic.MCONE_ID)
            fName = sprintf('AchromaticLMconeIsolatingSTFanalysis_RGC_%d_MconeCenter_%s.pdf', theRGCindex, strrep(opticsString, ', ', '_'));
        end

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));

        if (~isempty(pdfDirectory))
            NicePlot.exportFigToPDF(fullfile(pdfDirectory,fName), hFig, 300);
        end

    end % iRGC

    videoOBJ.close();

end

function [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight] = analyzeSurroundConeMixing(theMRGCmosaic, theRGCindex)

    [~, ~, theCenterMajorityConeType, ~] = theMRGCmosaic.centerConeTypeWeights(theRGCindex);

    % Retrieve this cell's  center cone indices and weights
    connectivityVector = full(squeeze(theMRGCmosaic.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    centerConeIndicesForCurrentRGC = find(connectivityVector > 0.0001);
    centerConeWeightsForCurrentRGC = connectivityVector(centerConeIndicesForCurrentRGC);

    
    % Retrieve this cell's surround cone indices and weights
    connectivityVector = full(squeeze(theMRGCmosaic.rgcRFsurroundConePoolingMatrix(:, theRGCindex)));
    surroundConeIndicesForCurrentRGC = find(connectivityVector > 0.0001);
    surroundConeWeightsForCurrentRGC = connectivityVector(surroundConeIndicesForCurrentRGC);

    idx = find(theMRGCmosaic.inputConeMosaic.coneTypes(centerConeIndicesForCurrentRGC) == cMosaic.LCONE_ID);
    netCenterLconeWeight = sum(centerConeWeightsForCurrentRGC(idx));

    idx = find(theMRGCmosaic.inputConeMosaic.coneTypes(centerConeIndicesForCurrentRGC) == cMosaic.MCONE_ID);
    netCenterMconeWeight = sum(centerConeWeightsForCurrentRGC(idx));

    idx = find(theMRGCmosaic.inputConeMosaic.coneTypes(surroundConeIndicesForCurrentRGC) == cMosaic.LCONE_ID);
    netSurroundLconeWeight = sum(surroundConeWeightsForCurrentRGC(idx));

    idx = find(theMRGCmosaic.inputConeMosaic.coneTypes(surroundConeIndicesForCurrentRGC) == cMosaic.MCONE_ID);
    netSurroundMconeWeight = sum(surroundConeWeightsForCurrentRGC(idx));
end

