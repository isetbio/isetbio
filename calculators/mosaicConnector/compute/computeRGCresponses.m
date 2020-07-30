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
    coVisualizedRetinalStimulusData = struct(...
        'spatialFrequency', p.Results.coVisualizedRetinalStimulusSpatialFrequency, ...
        'coneContrast', p.Results.coVisualizedRetinalStimulusConeContrast, ...
        'frameIndex', 1);
    
    % Assemble the null response filename
    [~, ~, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
    theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix, PolansSubjectID);
        
    % Open data file for read-only
    mFile = matfile(fullfile(saveDir,theNullResponseFileName), 'Writable', false);
    
    % Render retinal stimulus info
    coVisualizedRetinalStimulus = renderRetinalStimFigures(mFile, spatialFrequenciesCPD, ...
        coVisualizeRetinalStimulusWithMosaics, coVisualizedRetinalStimulusData, ...
        visualizeRetinalContrasts, visualizeMeanConeMosaicResponseAsAMovie, ...
        saveDir, runParams, LMScontrast, opticsPostFix, PolansSubjectID, ...
        figExportsDir);
    
    % Retrieve time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
    
    % Retrieve the null presynaptic responses
    switch presynapticSignal
        case 'isomerizations'
            theNullPresynapticResponses = mFile.isomerizationsNull;
        case   'photocurrents'
            theNullPresynapticResponses = mFile.photocurrentsNull;
        otherwise
            error('Unknown presynaptic signal: ''%s''.', presynapticSignal)
    end

    % Compute subregion responses to null stimulus
    [centerNullStimulusMeanResponse, surroundNullStimulusMeanResponse] = computeSubregionResponses(...
        theConeMosaic, theMidgetRGCmosaic.centerWeights, theMidgetRGCmosaic.surroundWeights, theNullPresynapticResponses);
       
    % Compute the RGC responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Assemble the test response filename
        theTestResponseFileName = sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID), gaborSpatialFrequencyCPD);
        
        % Open data file for read-only
        mFile = matfile(fullfile(saveDir, theTestResponseFileName), 'Writable', false);
       
        % Load the presynaptic input signals
        switch presynapticSignal
            case 'isomerizations'
                thePresynapticResponses = mFile.isomerizations;
            case   'photocurrents'
                thePresynapticResponses = mFile.photocurrents;
            otherwise
                error('UNknown presynaptic signal: ''%s''.', presynapticSignal)
        end
        
        % Compute the center and the surround responses
        fprintf('\nComputing RGC responses ...');
        tic
        
        % Compute center and surround subregion responses to the test stimulus
        [cRtest, sRtest] = computeSubregionResponses(theConeMosaic, theMidgetRGCmosaic.centerWeights, theMidgetRGCmosaic.surroundWeights, thePresynapticResponses);
            
        % Preallocate memory
        if (sfIndex == 1)
            instancesNum = size(cRtest,1);
            rgcsNum = size(cRtest,2);
            timeBins = size(cRtest,3);
            centerTestStimulusResponseInstances = zeros(numel(spatialFrequenciesCPD), instancesNum, rgcsNum, timeBins);
            surroundTestStimulusResponseInstances = centerTestStimulusResponseInstances;
        end
        
        centerTestStimulusResponseInstances(sfIndex,:,:,:) = cRtest;
        surroundTestStimulusResponseInstances(sfIndex,:,:,:) = sRtest;
        
        fprintf('Done in %2.1f minutes\n', toc/60);
    end % sfIndex
    
    
    % Center-surround integrated response instances
    integratedTestStimulusResponseInstances = centerTestStimulusResponseInstances - surroundTestStimulusResponseInstances;
    
    % Integrated response instance modulations (deviation from background response)
    integratedNullStimulusMeanResponse = centerNullStimulusMeanResponse - surroundNullStimulusMeanResponse;
    
    coneIsomerizationsDeltaPerSpikePerSecond = 600*1000.0;  % so many R*/sec lead to an RGC response of 1 spikes/sec above baseline
    responseTimeBin = responseTimeAxis(2)-responseTimeAxis(1);
    
    % Integrated responses in spikes/sec
    integratedResponseInstanceSpikesPerSec = isomerizationDeltasToSpikesPerSecond(...
            integratedTestStimulusResponseInstances, integratedNullStimulusMeanResponse, ...
            coneIsomerizationsDeltaPerSpikePerSecond, responseTimeBin);
    

    % Center responses in spikes/sec
    centerResponseInstances = isomerizationDeltasToSpikesPerSecond(...
            centerTestStimulusResponseInstances, centerNullStimulusMeanResponse, ...
            coneIsomerizationsDeltaPerSpikePerSecond, responseTimeBin);
        
    % Surround responses in spikes/sec
    surroundResponseInstances = isomerizationDeltasToSpikesPerSecond(...
            surroundTestStimulusResponseInstances, surroundNullStimulusMeanResponse, ...
            coneIsomerizationsDeltaPerSpikePerSecond, responseTimeBin);
        
    % Mean (over instances) integrated responses in spikes/sec
    integratedResponsesMean = squeeze(mean(integratedResponseInstanceSpikesPerSec,2));
    integratedResponsesStDev = squeeze(std(integratedResponseInstanceSpikesPerSec,0,2));

    
    % Report spikes per seconds
    for iRGC = 1:rgcsNum
        r = squeeze(integratedResponseInstanceSpikesPerSec(:, :, iRGC, :));
        r = squeeze(mean(r,2));
        fprintf('Max mean integrated response for RGC %d: %2.2f spikes/sec\n', iRGC, max(r(:)));
    end
    
    maxSpikeRateModulation = 200;
    maxSpikeRateModulationForComponents = 500;
    
    % Visualize the center and surround response components for the targeted RGCs
    if (visualizeResponseComponents)
        exportFig = true;
        for iTargetRGC = 1:numel(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves)
             visualizeResponseComponentsForTargetRGC(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC), ...
                 responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
                 spatialFrequenciesCPD, maxSpikeRateModulationForComponents, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir);  
        end
    end
    
    
    % Fit the DoG model to all SF tuning curves
    exportFig = true;
    [patchDogParams, responseAmplitude, spatialFrequenciesCPDHR, responseAmplitudeHR, responseTimeAxisHR, fittedResponsesHR] = ...
        computeAndFitDOGmodelToSpatialFrequencyTuningCurves(responseTimeAxis, integratedResponsesMean, integratedResponsesStDev,...
        maxSpikeRateModulation, stimSpatialParams, stimTemporalParams, spatialFrequenciesCPD, visualizeAllSpatialFrequencyTuningCurves, ...
        LMScontrast,  opticsPostFix, PolansSubjectID, exportFig, figExportsDir, ...
        'synthParams', theMidgetRGCmosaic.synthesizedRFParams.visual);
    
    % Generate figures showing temporal responses, SF tuning curves, and the DOG parameter statistics for the patch
    renderRGCanalysesFigures(patchDogParams, spatialFrequenciesCPDHR, responseAmplitudeHR, spatialFrequenciesCPD, responseAmplitude, ...
        responseTimeAxis, integratedResponsesMean, responseTimeAxisHR, fittedResponsesHR, ...
        maxSpikeRateModulation, theConeMosaic, runParams, theMidgetRGCmosaic, ...
        visualizePatchStatistics, visualizeRGCTemporalResponsesAtRGCPositions, visualizeRGCSFTuningsAtRGCPositions, ...
        targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves, coVisualizedRetinalStimulus, ...
        LMScontrast,  opticsPostFix, PolansSubjectID, exportFig, figExportsDir);
end

function responseInstanceModulationsSpikesPerSec = isomerizationDeltasToSpikesPerSecond(...
         responseInstances, nullResponses, ...
         coneIsomerizationsDeltaPerSpikePerSecond, responseTimeBinSeconds)
     
     sfsNum = size(responseInstances,1);
     instancesNum = size(responseInstances,2);
     rgcsNum = size(responseInstances,3);
     timeBins = size(responseInstances,4);
    
     % Subtract null responses (last time bin - useful for photocurrent response)
     nullResponses = reshape(nullResponses, [1 1 rgcsNum timeBins]);
     responseInstanceModulations = bsxfun(@minus, responseInstances, nullResponses(1,1,:,end));
     
     % Transform delta isomerizations/timebin to spikes/sec
     coneIsomerizationsDeltaPerTimeBin = coneIsomerizationsDeltaPerSpikePerSecond * responseTimeBinSeconds;
     responseInstanceModulationsSpikesPerSec = responseInstanceModulations / coneIsomerizationsDeltaPerTimeBin;
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
      
