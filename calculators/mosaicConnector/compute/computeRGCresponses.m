function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    presynapticSignal, spatialFrequenciesCPD, LMScontrast, stimSpatialParams, stimTemporalParams, ...
    saveDir, figExportsDir, ...
    visualizeRGCTemporalResponsesAtRGCPositions, visualizeRGCSFTuningsAtRGCPositions, ...
    visualizeAllSpatialFrequencyTuningCurves, visualizeResponseComponents, ...
    visualizeRetinalContrasts, visualizeMeanConeMosaicResponseAsAMovie, ...
    targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves, visualizePatchStatistics, ...
    varargin)
    
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Parse input
    p = inputParser;
    p.addParameter('coVisualizeRetinalStimulusWithMosaics', true, @islogical);
    p.addParameter('coVisualizedRetinalStimulusSpatialFrequency', 30, @isscalar);
    p.addParameter('coVisualizedRetinalStimulusConeContrast', LCONE_ID, @(x)(ismember(x, [LCONE_ID MCONE_ID SCONE_ID])));
    
    p.parse(varargin{:});
    coVisualizeRetinalStimulusWithMosaics = p.Results.coVisualizeRetinalStimulusWithMosaics;
    coVisualizedRetinalStimulusSpatialFrequency = p.Results.coVisualizedRetinalStimulusSpatialFrequency;
    coVisualizedRetinalStimulusConeContrast = p.Results.coVisualizedRetinalStimulusConeContrast;
    
    % Load the null presynaptic responses
    % Assemble the null response filename
    [~, ~, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
  
    theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix, PolansSubjectID);
        
    % Open data file for read-only
    mFile = matfile(fullfile(saveDir,theNullResponseFileName), 'Writable', false);
    switch presynapticSignal
        case 'isomerizations'
            theNullPresynapticResponses = mFile.isomerizationsNull;
        case   'photocurrents'
            theNullPresynapticResponses = mFile.photocurrentsNull;
        otherwise
            error('Unknown presynaptic signal: ''%s''.', presynapticSignal)
    end

    retinalStimulusDataIsAvailable = mFile.saveRetinalStimulusSequence;
    cornealStimulusDataIsAvailable = mFile.saveCornealStimulusSequence;
    
    % Retrieve background retinal LMS excitations
    if (coVisualizeRetinalStimulusWithMosaics && retinalStimulusDataIsAvailable)
        theRetinalLMSexcitationsSequenceNull = mFile.theRetinalLMSexcitationsSequence;
        backgroundLMSexcitations = mean(mean(mean(theRetinalLMSexcitationsSequenceNull,1),2),3);
    end
    
    % Retrieve time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
          
    % Stimulus to be co-visualized under the cone/RGC mosaic in target RGCs
    coVisualizedRetinalStimulus = [];
    
    % Compute the RGC responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Assemble the test response filename
        theTestResponseFileName = sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID), gaborSpatialFrequencyCPD);
        
        % Open data file for read-only
        mFile = matfile(fullfile(saveDir, theTestResponseFileName), 'Writable', false);
        
        if (coVisualizeRetinalStimulusWithMosaics && retinalStimulusDataIsAvailable)
            switch (coVisualizedRetinalStimulusConeContrast)
                case LCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 1;  % visualize the L-cone contrast
                case MCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 2;  % visualize the M-cone contrast
                case SCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 3;  % visualize the S-cone contrast
            end
            [~,visualizedRetinalSFindex] = min(abs(spatialFrequenciesCPD-coVisualizedRetinalStimulusSpatialFrequency));
            fprintf('Will superimpose retinal stimulus at %2.1f c/deg', spatialFrequenciesCPD(visualizedRetinalSFindex));
    
            if (sfIndex == visualizedRetinalSFindex)
                frameIndex = 5;      
                theRetinalLMSexcitations = mFile.theRetinalLMSexcitationsSequence(frameIndex,:,:,:);
                theRetinalLMScontrast = bsxfun(@times, bsxfun(@minus, theRetinalLMSexcitations, backgroundLMSexcitations), 1./backgroundLMSexcitations);
                coVisualizedRetinalStimulus.retinalContrastImage = squeeze(theRetinalLMScontrast(1,:,:,visualizedRetinalStimulusConeContrastIndex));
                coVisualizedRetinalStimulus.spatialSupportMicrons = mFile.retinalSpatialSupportMicrons;
            end
        end
        
        
        if (visualizeRetinalContrasts && retinalStimulusDataIsAvailable)
            % The test retinal LMS excitations 
            theRetinalLMSexcitationsSequence = mFile.theRetinalLMSexcitationsSequence;

            % Compute retinal contrast
            theRetinalLMScontrastSequence = bsxfun(@minus, theRetinalLMSexcitationsSequence, backgroundLMSexcitations);
            theRetinalLMScontrastSequence = bsxfun(@times, theRetinalLMScontrastSequence, 1./backgroundLMSexcitations);

            % Extract retinal L and M contrast over central region
            colsNum = size(theRetinalLMScontrastSequence,3);
            margin = round(colsNum/5);
            colsToUse = margin:2:(colsNum-margin);
            colsToUseNum = numel(colsToUse);
            midRow = round(0.5*size(theRetinalLMScontrastSequence,2));
            roiLcontrast(sfIndex,:) = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,1)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);
            roiMcontrast(sfIndex,:) = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,2)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);

            visualizeRetinalContrastProfiles(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast, figExportsDir);
        end

        
        % Load the presynaptic input signals
        switch presynapticSignal
            case 'isomerizations'
                thePresynapticResponses = mFile.isomerizations;
            case   'photocurrents'
                thePresynapticResponses = mFile.photocurrents;
            otherwise
                error('UNknown presynaptic signal: ''%s''.', presynapticSignal)
        end
        
        
        % Compute the differential response instances
        %thePresynapticResponses = bsxfun(@minus, thePresynapticResponses, theNullPresynapticResponses(:,:,end));
        
        % Compute the center and the surround responses
        fprintf('\nComputing RGC responses ...');
        tic
        
        % Compute subregion responses
        [cR, sR] = computeSubregionResponses(theConeMosaic, theMidgetRGCmosaic.centerWeights, theMidgetRGCmosaic.surroundWeights, thePresynapticResponses);
            
        if (sfIndex == 1)
            % Preallocate memory
            centerResponseInstances = zeros(numel(spatialFrequenciesCPD), size(cR,1), size(cR,2), size(cR,3));
            surroundResponseInstances = centerResponseInstances;
        end
        
        centerResponseInstances(sfIndex,:,:,:) = cR;
        surroundResponseInstances(sfIndex,:,:,:) = sR;
        
        fprintf('Done in %2.1f minutes\n', toc/60);
        
        % Display the mean presynaptic response
        if (visualizeMeanConeMosaicResponseAsAMovie && cornealStimulusDataIsAvailable)
            % Load the stimulus RGB sequence
            theStimulusRGBsequence = mFile.theStimulusSRGBsequence;
            theMeanPresynapticResponses = squeeze(mean(thePresynapticResponses,1));
            visualizeStimResponseMovie(responseTimeAxis, stimulusTimeAxis, stimSpatialParams, theStimulusRGBsequence, ...
                theMeanPresynapticResponses, theConeMosaic);
        end
    end % sfIndex
    
    if (visualizeRetinalContrasts)
        visualizeRetinalLMcontrastCorrelation(spatialFrequenciesCPD, roiLcontrast, roiMcontrast, LMScontrast, figExportsDir)
    end
    
    % Compute mean over all instances integrated response. This is used
    % to measure the F tuning below.
    integratedResponseInstances = centerResponseInstances-surroundResponseInstances;
   
    % Convert to spikes/sec
    maxSpikeRate = 200;
    m99 = prctile(integratedResponseInstances(:), 99);
   
    integratedResponseInstances = integratedResponseInstances / m99 * maxSpikeRate;
    integratedResponsesMean = squeeze(mean(integratedResponseInstances,2));
    integratedResponsesStDev = squeeze(std(integratedResponseInstances,0,2));
        
    mTotal = max([max(abs(centerResponseInstances(:))) max(abs(surroundResponseInstances(:)))]);
    centerResponseInstances = centerResponseInstances / mTotal * maxSpikeRate;
    surroundResponseInstances = surroundResponseInstances / mTotal * maxSpikeRate;
   
    labelCells = true;
    if (visualizeResponseComponents)
        for iTargetRGC = 1:numel(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves)
             visualizeResponseComponentsForTargetRGC(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC), responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
                spatialFrequenciesCPD, maxSpikeRate, LMScontrast, figExportsDir);  
        end
    end
    
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        
        % Compute response tuning
        switch (stimSpatialParams.type)
            case 'driftingGrating'
                visualizeIndividualFits = ~true;
                exportFig = ~true;
                % Compute modulation of the response at the fundamental temporal frequency
                [responseAmplitude(:, sfIndex), responseAmplitudeSE(:, sfIndex), responsePhase(:, sfIndex), ...
                    responseTimeAxisHR, fittedResponsesHR(sfIndex,:,:)] = fitSinusoidalModulationToResponseTimeCourse(...
                    squeeze(integratedResponsesMean(sfIndex,:,:)), ...
                    squeeze(integratedResponsesStDev(sfIndex,:,:)), ...
                    responseTimeAxis, ...
                    stimTemporalParams.temporalFrequencyHz, ...
                    spatialFrequenciesCPD(sfIndex), maxSpikeRate, ...
                    visualizeIndividualFits, LMScontrast, [], exportFig, figExportsDir);
            otherwise
                error('Unknown stimulus type: ''%''.', stimulusSpatialParams.type)
        end
        
    end % sfIndex
    
    % Fit the DoG model to the cell SF tuning
    maxSpikeRatModulation = 26;
    initialParams = [];
    visualizeIndividualFits = false; exportFig = true;
    [~,~,~, meanParams] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, responseAmplitudeSE, initialParams, ...
        maxSpikeRatModulation, LMScontrast, visualizeIndividualFits, exportFig, '');
    
    % Second fit
    visualizeIndividualFits = visualizeAllSpatialFrequencyTuningCurves; exportFig = true;
    initialParams = meanParams;
    [patchDogParams,spatialFrequenciesCPDHR, responseAmplitudeHR] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, responseAmplitudeSE, initialParams, ...
        maxSpikeRatModulation, LMScontrast, visualizeIndividualFits, exportFig, figExportsDir);

    if (visualizePatchStatistics)
        % Visualize data to contrast with Cronner and Kaplan data
        RGCpositionsMicrons = determineRGCPositionsFromCenterInputs(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, theMidgetRGCmosaic.centerWeights);
        RGCeccentricityDegs = WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCpositionsMicrons.^2,2))/1000.0);
        visualizePatchStatsDerivedFromSFcurves(patchDogParams, RGCeccentricityDegs);
    end
    
    %   Visualize the temporal response of each RGC at the RGC's location
    if (visualizeRGCTemporalResponsesAtRGCPositions)
        for sfIndex = 0:numel(spatialFrequenciesCPD)
            exportFig = true;
            superimposedRetinalStimulus = [];
            plotXaxisScaling = 'linear';
            plotType = 'TimeResponse';
            figureName = sprintf('%2.1fcpdResponse', spatialFrequenciesCPD(sfIndex));
            visualizeRGCmosaicWithResponses(100+sfIndex, theConeMosaic, plotXaxisScaling, plotType, ...
               responseTimeAxis, squeeze(integratedResponsesMean(sfIndex,:,:)), ...
               responseTimeAxisHR, squeeze(fittedResponsesHR(sfIndex,:,:)), ...
               runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
               theMidgetRGCmosaic, 'centers', maxSpikeRate, ...
               superimposedRetinalStimulus, ....
               figureName, LMScontrast, opticsPostFix, PolansSubjectID, ...
               [], labelCells, ...
               exportFig, figExportsDir);
        end
    end
    
    if (visualizeRGCSFTuningsAtRGCPositions)
        % Visualize the response tuning of each RGC at the RGC's location
        exportFig = true;
        superimposedRetinalStimulus = [];
        plotXaxisScaling = 'log';
        plotType = 'SFtuning';
        figureName = 'SFtuningAll';
        visualizeRGCmosaicWithResponses(1000, theConeMosaic, plotXaxisScaling, plotType, ...
                    spatialFrequenciesCPD, responseAmplitude, ...
                    spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                    runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                    theMidgetRGCmosaic, 'centers', maxSpikeRatModulation, ...
                    superimposedRetinalStimulus, ....
                    figureName, LMScontrast, opticsPostFix, PolansSubjectID, ...
                    [], labelCells, ...
                    exportFig, figExportsDir);
    end
    
    if (~isempty(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves))
        for iTargetRGC = 1:numel(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves)
            figureName = sprintf('SFtuning%d', targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC));

            visualizeRGCmosaicWithResponses(1000+targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC), theConeMosaic, plotXaxisScaling, plotType, ...
                    spatialFrequenciesCPD, responseAmplitude, ...
                    spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                    runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                    theMidgetRGCmosaic,  'centers', maxSpikeRatModulation, ...
                    coVisualizedRetinalStimulus, ...
                    figureName, LMScontrast, opticsPostFix, PolansSubjectID, ...
                    targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC), false, ...
                    exportFig, figExportsDir);
        end
    end
    
    % Visualize the mosaics
    visualizeMosaics = ~true;
    if (visualizeMosaics)
        fprintf('\nVisualizing the RGC mosaic with the optical image ...');
        theOISequence = mFile.theOIsequence;
        theFirstOI = theOISequence.frameAtIndex(1);
        zLevels = [0.3 1];
       
        hFig = visualizeConeAndRGCmosaicsWithRetinalImage(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
            theMidgetRGCmosaic, zLevels, 'centers', theFirstOI); 
        plotlabOBJ.exportFig(hFig, 'pdf', sprintf('%s.mat',coneResponsesFileName), pwd());
        fprintf('Done !\n');
    end
end

function [responsesC, responsesS] = computeSubregionResponses(theConeMosaic, weightsC, weightsS, presynapticResponses)
    % Get dimensionalities
    [instancesNum, conesNum, timeBins] = size(presynapticResponses);
    rgcsNum = size(weightsC,2);
    
    % Form response matrix
    responsesC = zeros(instancesNum, rgcsNum, timeBins);
    responsesS = zeros(instancesNum, rgcsNum, timeBins);
    
    t = (0:(timeBins-1)) * theConeMosaic.integrationTime;
    tauC = 2.5/1000;
    centerIR = exp(-t/tauC);
    centerIR = centerIR/sum(centerIR);
    tauC = 25/1000;
    surroundIR = exp(-t/tauC);
    surroundIR = surroundIR/sum(surroundIR);
    
    for instanceIndex = 1:instancesNum
        % All presynaptic spatiotemporal responses for this instance
        instancePresynapticResponse = squeeze(presynapticResponses(instanceIndex,:,:));
        for iRGC = 1:rgcsNum
            % The RGC's weights
            iRGCweightsC = (full(squeeze(weightsC(:,iRGC))))';
            iRGCweightsS = (full(squeeze(weightsS(:,iRGC))))';
            
            % The RGC temporal response
            centerR = iRGCweightsC * instancePresynapticResponse;
            surroundR = iRGCweightsS * instancePresynapticResponse;
%             
%             centerR = conv(centerR, centerIR);
%             centerR = centerR(1:timeBins);
%             
%             surroundR = conv(surroundR, surroundIR);
%             surroundR = surroundR(1:timeBins);
            
            responsesC(instanceIndex,iRGC,:) = centerR;
            responsesS(instanceIndex,iRGC,:) = surroundR;
        end % iRGC
        
    end % instanceIndex
end
       
