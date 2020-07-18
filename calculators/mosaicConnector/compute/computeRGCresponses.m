function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    presynapticSignal, spatialFrequenciesCPD, LMScontrast, stimSpatialParams, stimTemporalParams, ...
    targetRGCs, saveDir, figExportsDir, visualizeAllSpatialFrequencyTuningCurves, visualizeResponseComponents)
    
    % Load the null presynaptic responses
    % Assemble the null response filename
    [~, ~, opticsPostFix] = mosaicsAndOpticsFileName(runParams);
    theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix);
        
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

    % Retrieve background LMS excitations
    theRetinalLMSexcitationsSequenceNull = mFile.theRetinalLMSexcitationsSequence;
    backgroundLMSexcitations = mean(mean(mean(theRetinalLMSexcitationsSequenceNull,1),2),3);
    
    % Retrieve time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
          
    visualizedRetinalStimulus = [];
    visualizedRetinalStimulusSpatialFrequencyCPD = 60;
    visualizedRetinalStimulusConeContrastIndex = 1;  % visualize the L-cone contrast
    [~,visualizedRetinalSFindex] = min(abs(spatialFrequenciesCPD-visualizedRetinalStimulusSpatialFrequencyCPD));
    
    fprintf('Will superimpose retinal stimulus at %2.1f c/deg', spatialFrequenciesCPD(visualizedRetinalSFindex));
    
    % Compute the RGC responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Assemble the test response filename
        theTestResponseFileName = sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast, opticsPostFix), gaborSpatialFrequencyCPD);
        
        % Open data file for read-only
        mFile = matfile(fullfile(saveDir, theTestResponseFileName), 'Writable', false);
        
        if (sfIndex == visualizedRetinalSFindex)
            frameIndex = 5;      
            size(mFile.theRetinalLMSexcitationsSequence)
            theRetinalLMSexcitations = mFile.theRetinalLMSexcitationsSequence(frameIndex,:,:,:);
            theRetinalLMScontrast = bsxfun(@times, bsxfun(@minus, theRetinalLMSexcitations, backgroundLMSexcitations), 1./backgroundLMSexcitations);
            visualizedRetinalStimulus.retinalContrastImage = squeeze(theRetinalLMScontrast(1,:,:,visualizedRetinalStimulusConeContrastIndex));
            visualizedRetinalStimulus.spatialSupportMicrons = mFile.retinalSpatialSupportMicrons;
        end

        visualizeRetinalContrasts = ~true;
        if (visualizeRetinalContrasts)     
            % The test retinal LMS excitations 
            %theOISRGBsequence = mFile.theOISRGBsequence;
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

            visualizeRetinalContrastSequence(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast, figExportsDir);
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
        showMeanMosaicResponseAsMovie = ~true;
        if (showMeanMosaicResponseAsMovie)
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
        for iTargetRGC = 1:numel(targetRGCs)
             visualizeResponseComponentsForTargetRGC(targetRGCs(iTargetRGC), responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
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
    
    % Fit the responses
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

    % Visualize data to contrast with Cronner and Kaplan data
    RGCpositionsMicrons = determineRGCPositionsFromCenterInputs(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, theMidgetRGCmosaic.centerWeights);
    RGCeccentricityDegs = WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCpositionsMicrons.^2,2))/1000.0);
    visualizePatchStatsDerivedFromSFcurves(patchDogParams, RGCeccentricityDegs);
    
%   Visualize the temporal response of each RGC at the RGC's location
    for sfIndex = 0:-1 %1:numel(spatialFrequenciesCPD)
        
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
           figureName, LMScontrast, [], labelCells, ...
           exportFig, figExportsDir);
    end
    
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
                figureName, LMScontrast, [], labelCells, ...
                exportFig, figExportsDir);
            
    for iTargetRGC = 1:numel(targetRGCs)
        figureName = sprintf('SFtuning%d', targetRGCs(iTargetRGC));

        visualizeRGCmosaicWithResponses(1000+targetRGCs(iTargetRGC), theConeMosaic, plotXaxisScaling, plotType, ...
                spatialFrequenciesCPD, responseAmplitude, ...
                spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                theMidgetRGCmosaic,  'centers', maxSpikeRatModulation, ...
                visualizedRetinalStimulus, ...
                figureName, LMScontrast, targetRGCs(iTargetRGC), false, ...
                exportFig, figExportsDir);
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
       
