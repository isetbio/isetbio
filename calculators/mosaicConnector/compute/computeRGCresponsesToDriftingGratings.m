function computeRGCresponsesToDriftingGratings(runParams, theConeMosaic, theConeMosaicMetaData, theMidgetRGCmosaic, ...
    presynapticSignal, LMScontrast, stimSpatialParams, stimTemporalParams, ...
    recomputeEccBasedPhotocurrentResponses, ...
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
    p.addParameter('coVisualizedRetinalStimulusSpatialFrequencyCPD', 30, @isscalar);
    p.addParameter('coVisualizedRetinalStimulusConeContrast', LCONE_ID, @(x)(ismember(x, [LCONE_ID MCONE_ID SCONE_ID])));
    
    p.parse(varargin{:});
    coVisualizeRetinalStimulusWithMosaics = p.Results.coVisualizeRetinalStimulusWithMosaics;
    coVisualizedRetinalStimulusData = struct(...
        'spatialFrequency', p.Results.coVisualizedRetinalStimulusSpatialFrequencyCPD, ...
        'coneContrast', p.Results.coVisualizedRetinalStimulusConeContrast, ...
        'frameIndex', 1);
    
    % Assemble the null response filename
    [~, ~, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
    theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix, PolansSubjectID);
        

    % Open data file for read-only
    mFile = matfile(fullfile(saveDir,theNullResponseFileName), 'Writable', false);
    
    % Extract tested spatial frequencies
    spatialFrequenciesCPD = stimSpatialParams.testedSpatialFrequenciesCPD;
     
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
    theIsomerizationsNull = mFile.isomerizationsNull;
    
    switch presynapticSignal
        case 'isomerizations'
            theNullPresynapticResponses = mFile.isomerizationsNull;
        case   'photocurrents'
            if (recomputeEccBasedPhotocurrentResponses)
                theNullPresynapticResponses = computeEccBasedPhotocurrentsFromIsomerizations(...
                        theConeMosaic, theConeMosaicMetaData, mFile.isomerizationsNull, theIsomerizationsNull, runParams.rgcMosaicPatchEccMicrons);
            else
                theNullPresynapticResponses = mFile.photocurrentsNull;
            end
            
        otherwise
            error('Unknown presynaptic signal: ''%s''.', presynapticSignal)
    end
    
    
    
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
                    % Using Weber response representation ensures that the cone
                    % responses are proportional to cone Weber contrast. 
                    % That way a 10% L+M cone contrast will active L-cone center RGCs
                    % to the same extent as it will activate M-cone center RGCs. If we
                    % do not use this the relative response of L-cone center to M-cone center RGCs
                    % depends on the mean cone response to the background.
                    % Phototransductios does this transformation so that
                    % pulses of same contrast have the same photocurrent modulation
                    % independent on the background activation. But since we are
                    % working on the cone isomerizations signal here, we have to do
                    % this contrast transformation to have RGC activations that are
                    % proportional to stimuls contrast.
                useWeberResponseRepresentation = true;
                thePresynapticResponses = mFile.isomerizations;
                
            case 'photocurrents'
                useWeberResponseRepresentation = ~true;
                if (recomputeEccBasedPhotocurrentResponses)
                    fprintf('\nComputing photocurrent responses ...');
                    tic
                    thePresynapticResponses = computeEccBasedPhotocurrentsFromIsomerizations(...
                        theConeMosaic, theConeMosaicMetaData, mFile.isomerizations, theIsomerizationsNull, runParams.rgcMosaicPatchEccMicrons);
                    fprintf('Done in %2.1f minutes\n', toc/60);
                else
                    thePresynapticResponses = mFile.photocurrents;
                    % Subtract the steady state background photocurrents
                    thePresynapticResponses = bsxfun(@minus, thePresynapticResponses, theNullPresynapticResponses(1,:,end));
                end
            otherwise
                error('UNknown presynaptic signal: ''%s''.', presynapticSignal)
        end
           
        
        % Compute the center and the surround responses
        fprintf('\nComputing RGC responses ...');
        tic
        if (useWeberResponseRepresentation)
             % Compute the pre-spatial responses as Weber responses (test-null)/null
             preSpatialIntegrationResponses = bsxfun(@minus, thePresynapticResponses, theNullPresynapticResponses);
             preSpatialIntegrationResponses = bsxfun(@times, preSpatialIntegrationResponses, 1./(theNullPresynapticResponses+eps));
             preSpatialIntergrationConversionFactorToSpikePerSecond = 5;
        else
            % Compute the pre-spatial integration responses as differential responses (test-null)
            preSpatialIntegrationResponses = thePresynapticResponses;
            % Conversion to spikes/sec
            preSpatialIntergrationConversionFactorToSpikePerSecond = 10.0;
        end
        
        % Spatial intergration using center/surround weights
        [cRtest, sRtest] = computeSubregionResponses(theConeMosaic, theMidgetRGCmosaic.centerWeights, ...
            theMidgetRGCmosaic.surroundWeights, preSpatialIntegrationResponses);
            
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
    integratedTestStimulusResponseInstances = centerTestStimulusResponseInstances + surroundTestStimulusResponseInstances;
    
    % Integrated responses in spikes/sec
    integratedResponseInstanceSpikesPerSec = ...
        integratedTestStimulusResponseInstances / preSpatialIntergrationConversionFactorToSpikePerSecond;
    
    % Center responses in spikes/sec
    centerResponseInstances = ...
        centerTestStimulusResponseInstances / preSpatialIntergrationConversionFactorToSpikePerSecond;
        
    % Surround responses in spikes/sec
    surroundResponseInstances = ...
        surroundTestStimulusResponseInstances/ preSpatialIntergrationConversionFactorToSpikePerSecond;
        
    % Mean (over instances) integrated responses in spikes/sec
    integratedResponsesMean = squeeze(mean(integratedResponseInstanceSpikesPerSec,2));
    integratedResponsesStDev = squeeze(std(integratedResponseInstanceSpikesPerSec,0,2));

    % Report spikes per second for each RGC
    maxResponse = zeros(1, rgcsNum);
    for iRGC = 1:rgcsNum
        r = squeeze(integratedResponseInstanceSpikesPerSec(:, :, iRGC, :));
        r = squeeze(mean(r,2));
        maxResponse(iRGC) = max(r(:));
        fprintf('Max mean integrated response for RGC %d: %2.2f spikes/sec\n', iRGC, maxResponse(iRGC));
    end
    
    % Select SpikeRate max and ticks
    maxAllResponses = prctile(integratedResponsesMean(:), 99);
   
    if (maxAllResponses > 300)
        maxSpikeRateModulation = ceil(maxAllResponses/50)*50;
        spikeRateTicks = 0:50:maxSpikeRateModulation;
    elseif (maxAllResponses > 200)
        maxSpikeRateModulation = 300;
        spikeRateTicks = 0:50:maxSpikeRateModulation;
    elseif (maxAllResponses > 100)
        maxSpikeRateModulation = 200;
        spikeRateTicks = 0:40:maxSpikeRateModulation;
    elseif (maxAllResponses > 50)
        maxSpikeRateModulation = 100;
        spikeRateTicks = 0:20:maxSpikeRateModulation;
    elseif (maxAllResponses > 25)
        maxSpikeRateModulation = 50;
        spikeRateTicks = 0:10:maxSpikeRateModulation;
    else
        maxSpikeRateModulation = 25;
        spikeRateTicks = 0:5:maxSpikeRateModulation;
    end

    % Visualize the center and surround response components for the targeted RGCs
    if (visualizeResponseComponents)
        exportFig = true;
        for iTargetRGC = 1:numel(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves)
             visualizeResponseComponentsForTargetRGC(targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves(iTargetRGC), ...
                 responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
                 spatialFrequenciesCPD, maxSpikeRateModulation, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir);  
        end
    end
    
    % Fit the DoG model to all SF tuning curves
    exportFig = true;
    [patchDogParams, modelFitted, responseAmplitude, spatialFrequenciesCPDHR, ...
     responseAmplitudeHR, responseTimeAxisHR, fittedResponsesHR] = fitSpatialTransferFunctionData(responseTimeAxis, ...
        integratedResponsesMean, integratedResponsesStDev,...
        maxSpikeRateModulation, spikeRateTicks, stimSpatialParams, stimTemporalParams, LMScontrast, ...
        spatialFrequenciesCPD, visualizeAllSpatialFrequencyTuningCurves, ...
        targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves, ...
        opticsPostFix, PolansSubjectID, exportFig, figExportsDir, ...
        'synthParams', theMidgetRGCmosaic.synthesizedRFParams.visual);
    
    % Export the computed spatial transfer functions
    theSpatialTransferFunctionsFileName = fullfile(saveDir, spatialTransferFunctionsFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID));
    if (LMScontrast(1) == LMScontrast(2))
         save(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude', ...
             'theConeMosaic', 'theMidgetRGCmosaic');
    else
        save(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude');
    end
    
    fprintf('Spatial transfer functions for  [%2.1f, %2.1f %2.1f] stimulus saved to %s.\n', ...
        LMScontrast(1), LMScontrast(2), LMScontrast(3), theSpatialTransferFunctionsFileName);

    % Generate figures showing temporal responses, SF tuning curves, and the DOG parameter statistics for the patch
    renderRGCanalysesFigures(patchDogParams, modelFitted, ...
        spatialFrequenciesCPDHR, responseAmplitudeHR, ...
        spatialFrequenciesCPD, responseAmplitude, ...
        responseTimeAxis, integratedResponsesMean, responseTimeAxisHR, fittedResponsesHR, ...
        maxSpikeRateModulation, spikeRateTicks, theConeMosaic, runParams, theMidgetRGCmosaic, ...
        visualizePatchStatistics, visualizeRGCTemporalResponsesAtRGCPositions, visualizeRGCSFTuningsAtRGCPositions, ...
        targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves, coVisualizedRetinalStimulus, ...
        LMScontrast,  opticsPostFix, PolansSubjectID, exportFig, figExportsDir);
end


function [responsesC, responsesS] = computeSubregionResponses(theConeMosaic, weightsC, weightsS, preSpatialIntegrationResponses)
    % Get dimensionalities
    [instancesNum, conesNum, timeBins] = size(preSpatialIntegrationResponses);
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
        instancePreSpatialIntegrationResponse = squeeze(preSpatialIntegrationResponses(instanceIndex,:,:));
        for iRGC = 1:rgcsNum
            % The RGC's weights
            iRGCweightsC = (full(squeeze(weightsC(:,iRGC))))';
            iRGCweightsS = (full(squeeze(weightsS(:,iRGC))))';
            
            % The RGC temporal response
            centerR = iRGCweightsC * instancePreSpatialIntegrationResponse;
            surroundR = iRGCweightsS * instancePreSpatialIntegrationResponse;
%             
%             centerR = conv(centerR, centerIR);
%             centerR = centerR(1:timeBins);
%             
%             surroundR = conv(surroundR, surroundIR);
%             surroundR = surroundR(1:timeBins);
            
            responsesC(instanceIndex,iRGC,:) = centerR;
            responsesS(instanceIndex,iRGC,:) = -surroundR;
        end % iRGC
        
    end % instanceIndex
end
      
