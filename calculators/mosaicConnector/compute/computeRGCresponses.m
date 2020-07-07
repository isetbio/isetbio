function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    presynapticSignal, spatialFrequenciesCPD, LMScontrast, stimSpatialParams, stimTemporalParams, saveDir)
    
    % Load the null presynaptic responses
    mFile = matfile(fullfile(saveDir,nullResponseFilename(runParams)), 'Writable', false);
    switch presynapticSignal
        case 'isomerizations'
            theNullPresynapticResponses = mFile.isomerizationsNull;
        case   'photocurrents'
            theNullPresynapticResponses = mFile.photocurrentsNull;
        otherwise
            error('UNknown presynaptic signal: ''%s''.', presynapticSignal)
    end

    % Retrieve time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
           
    % Compute the RGC responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        dataFile = fullfile(saveDir, sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast), gaborSpatialFrequencyCPD));
        mFile = matfile(dataFile, 'Writable', false);
    
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
        thePresynapticResponses = bsxfun(@minus, thePresynapticResponses, theNullPresynapticResponses(:,:,end));
        
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
        showMeanMosaicResponseAsMovie = ~true;
        if (showMeanMosaicResponseAsMovie)
            % Load the stimulus RGB sequence
            theStimulusRGBsequence = mFile.theStimulusRGBsequence;
            theMeanPresynapticResponses = squeeze(mean(thePresynapticResponses,1));
            visualizeStimResponseMovie(responseTimeAxis, stimulusTimeAxis, theStimulusRGBsequence, theMeanPresynapticResponses, theConeMosaic);
        end
    end % sfIndex
    
    
   % Compute mean over all instances integrated response. This is used
   % to measure the F tuning below.
   integratedResponseInstances = centerResponseInstances-surroundResponseInstances;
   
   % Convert to spikes/sec
   maxSpikeRate = 200;
   m99 = prctile(integratedResponseInstances(:), 99);
   
   integratedResponseInstances = integratedResponseInstances / m99 * maxSpikeRate;
   integratedResponsesMean = squeeze(mean(integratedResponseInstances,2));
   integratedResponsesStDev = squeeze(std(integratedResponseInstances,0,2));
        
   for sfIndex = 1:numel(spatialFrequenciesCPD)
        
        % Compute response tuning
        switch (stimSpatialParams.type)
            case 'driftingGrating'
                visualizeFits = ~true;
                % Compute modulation of the response at the fundamental temporal frequency
                [responseAmplitude(:, sfIndex), responsePhase(:, sfIndex), ...
                    responseTimeAxisHR, fittedResponsesHR(sfIndex,:,:)] = fitSinusoidalModulationToResponseTimeCourse(...
                    squeeze(integratedResponsesMean(sfIndex,:,:)), ...
                    squeeze(integratedResponsesStDev(sfIndex,:,:)), ...
                    responseTimeAxis, ...
                    stimTemporalParams.temporalFrequencyHz, ...
                    visualizeFits, spatialFrequenciesCPD(sfIndex), maxSpikeRate);
            otherwise
                error('Unknown stimulus type: ''%''.', stimulusSpatialParams.type)
        end
    end
    
    % Fit the responses
    plotFitResults = ~true; initialParams = [];
    [~,~,~, meanParams] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, plotFitResults, initialParams, maxSpikeRate);
    
    % Second fit
    plotFitResults = true; initialParams = meanParams;
    [patchDogParams,spatialFrequenciesCPDHR, responseAmplitudeHR, meanParams] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, plotFitResults, initialParams, maxSpikeRate);

    % Visualize data to contrast with Cronner and Kaplan data
    RGCpositionsMicrons = determineRGCPositionsFromCenterInputs(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, theMidgetRGCmosaic.centerWeights);
    RGCeccentricityDegs = WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCpositionsMicrons.^2,2))/1000.0);
    visualizePatchStatsDerivedFromSFcurves(patchDogParams, RGCeccentricityDegs);
    
    iRGC = 100;
    figNo = 123;
    visualizeResponseComponentsForSingleRGC(figNo, iRGC, responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
        responseTimeAxisHR, squeeze(fittedResponsesHR(:,iRGC,:)), spatialFrequenciesCPD, maxSpikeRate);  
    
    % Visualize the temporal response of each RGC at the RGC's location
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        zLevels = [0.3 1];
        visualizeRGCmosaicWithResponses(100+sfIndex, theConeMosaic, 'linear', ...
           responseTimeAxis, squeeze(integratedResponsesMean(sfIndex,:,:)), ...
           responseTimeAxisHR, squeeze(fittedResponsesHR(sfIndex,:,:)), ...
           runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
           theMidgetRGCmosaic, zLevels, 'centers', maxSpikeRate);    
    end
    
    % Visualize the response tuning of each RGC at the RGC's location
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, 'log', ...
                spatialFrequenciesCPD, responseAmplitude, ...
                spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                theMidgetRGCmosaic, zLevels, 'centers', maxSpikeRate); 
   
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
        
function  visualizeResponseComponentsForSingleRGC(figNo, iRGC, responseTimeAxis, centerResponses, surroundResponses, ...
    responseTimeAxisHR, fittedResponses, spatialFrequenciesCPD, maxSpikeRate)

    centerResponsesMean = squeeze(mean(centerResponses,2));
    surroundResponsesMean = squeeze(mean(surroundResponses,2));
    m = maxSpikeRate;

    figure(figNo); clf;
    for sfIndex = 1:numel(spatialFrequenciesCPD)    
        subplot(3,5,sfIndex);
        centerResponses = squeeze(centerResponsesMean(sfIndex,iRGC,:));
        plot(responseTimeAxis, centerResponses, 'rs-'); hold on;
        surroundResponses = squeeze(surroundResponsesMean(sfIndex,iRGC,:));
        plot(responseTimeAxis, surroundResponses, 'bs-');
        plot(responseTimeAxisHR, squeeze(fittedResponses(sfIndex,:)), 'k-', 'LineWidth', 1.5);
        set(gca, 'YLim', [-m m]);
    end
end
