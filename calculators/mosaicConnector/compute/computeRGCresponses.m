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
            error('Unknown presynaptic signal: ''%s''.', presynapticSignal)
    end

    
    visualizeRetinalContrasts = ~true;
    if (visualizeRetinalContrasts)
        % Retrieve the background retinal LMS excitations
        %theOISRGBsequence = mFile.theOISRGBsequence;
        theRetinalLMSexcitationsSequenceNull = mFile.theRetinalLMSexcitationsSequence;
        backgroundLMSexcitations = mean(mean(mean(theRetinalLMSexcitationsSequenceNull,1),2),3);
    end
    
    % Retrieve time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
           
    % Compute the RGC responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        dataFile = fullfile(saveDir, sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast), gaborSpatialFrequencyCPD));
        mFile = matfile(dataFile, 'Writable', false);
        
        
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
            colsToUse = margin:2:colsNum-margin;
            colsToUseNum = numel(colsToUse);
            midRow = round(0.5*size(theRetinalLMScontrastSequence,2));
            roiLcontrast(sfIndex,:) = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,1)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);
            roiMcontrast(sfIndex,:) = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,2)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);

            visualizeRetinalContrastSequence(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast);
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
        visualizeRetinalLMcontrastCorrelation(spatialFrequenciesCPD, roiLcontrast, roiMcontrast, LMScontrast)
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
    targetRGC = 1;  % 36 (M-only center, diagonal) 1 (M-only center, 2 vertical cones) 120 (L-only center, artifact)
    visualizeResponseComponentsForSingleRGC(targetRGC, responseTimeAxis, centerResponseInstances, surroundResponseInstances, ...
        spatialFrequenciesCPD, maxSpikeRate, LMScontrast);  
  
    
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
                    visualizeIndividualFits, exportFig, LMScontrast, targetRGC);
            otherwise
                error('Unknown stimulus type: ''%''.', stimulusSpatialParams.type)
        end
    end
    
    % Fit the responses
    maxSpikeRatModulation = 26;
    initialParams = [];
    [~,~,~, meanParams] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, responseAmplitudeSE, initialParams, ...
        maxSpikeRatModulation, false, false, LMScontrast);
    
    % Second fit
    visualizeIndividualFits = true; exportFig = true;
    initialParams = meanParams;
    [patchDogParams,spatialFrequenciesCPDHR, responseAmplitudeHR] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, responseAmplitudeSE, initialParams, ...
        maxSpikeRatModulation, visualizeIndividualFits, exportFig, LMScontrast);

    % Visualize data to contrast with Cronner and Kaplan data
    RGCpositionsMicrons = determineRGCPositionsFromCenterInputs(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, theMidgetRGCmosaic.centerWeights);
    RGCeccentricityDegs = WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCpositionsMicrons.^2,2))/1000.0);
    visualizePatchStatsDerivedFromSFcurves(patchDogParams, RGCeccentricityDegs);
    
    zLevels = [0.3 1];
     
   
%     % Visualize the temporal response of each RGC at the RGC's location
    for sfIndex = 0:-1 %1:numel(spatialFrequenciesCPD)
        
        exportFig = true;
        visualizeRGCmosaicWithResponses(100+sfIndex, theConeMosaic, 'linear', ...
           responseTimeAxis, squeeze(integratedResponsesMean(sfIndex,:,:)), ...
           responseTimeAxisHR, squeeze(fittedResponsesHR(sfIndex,:,:)), ...
           runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
           theMidgetRGCmosaic, zLevels, 'centers', maxSpikeRate, ...
           exportFig, sprintf('%2.1fcpdResponse', spatialFrequenciesCPD(sfIndex)), LMScontrast, [], labelCells);
    end
    
    % Visualize the response tuning of each RGC at the RGC's location
    exportFig = true;
    
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, 'log', ...
                spatialFrequenciesCPD, responseAmplitude, ...
                spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                theMidgetRGCmosaic, zLevels, 'centers', maxSpikeRatModulation, ...
                exportFig, 'SFtuningAll', LMScontrast, [], labelCells);
            
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, 'log', ...
                spatialFrequenciesCPD, responseAmplitude, ...
                spatialFrequenciesCPDHR, responseAmplitudeHR, ...
                runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                theMidgetRGCmosaic, zLevels, 'centers', maxSpikeRatModulation, ...
                exportFig, sprintf('SFtuning%d', targetRGC), LMScontrast, targetRGC, false);
   
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
    
    for iRGC = 1:rgcsNum
            % The RGC's weights
            iRGCweightsC = full(squeeze(weightsC(:,iRGC)));
            iRGCweightsS = full(squeeze(weightsS(:,iRGC)));
            iRatio(iRGC) = sum(iRGCweightsS(:))/sum(iRGCweightsC(:));
            fprintf(2,'RGC %d integrated S/C ratio from weights: %2.2f\n', iRGC, iRatio(iRGC));
    end
    fprintf(2, 'mean  integrated S/C ratio from weights: %2.2f\n', mean(iRatio));
    
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
        
function  visualizeResponseComponentsForSingleRGC(targetRGC, responseTimeAxis, centerResponses, surroundResponses, ...
    spatialFrequenciesCPD, maxSpikeRate, LMScontrast)

    plotlabOBJ = setupPlotLab(0, 16, 12);
    hFig = figure(123); clf;
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 3, ...
        'colsNum', 4, ...
        'leftMargin', 0.06, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.03, ...
        'bottomMargin', 0.05, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    centerResponsesMean = squeeze(mean(centerResponses,2));
    surroundResponsesMean = squeeze(mean(surroundResponses,2));
    
    for sfIndex = 1:numel(spatialFrequenciesCPD)  
        row = floor((sfIndex-1)/4)+1;
        col = mod(sfIndex-1,4)+1;
        ax = theAxesGrid{row,col};
        centerResponses = squeeze(centerResponsesMean(sfIndex,targetRGC,:));
        surroundResponses = squeeze(surroundResponsesMean(sfIndex,targetRGC,:));
        line(ax, responseTimeAxis, centerResponses, 'Color', [1 0 0], 'LineWidth', 1.5); hold on;
        line(ax, responseTimeAxis, surroundResponses, 'Color', [0 0 1], 'LineWidth', 1.5);
        set(ax, 'XLim', [0 0.5], 'YLim', maxSpikeRate*[0 1], 'XTick', 0:0.1:0.5, 'YTick',(0:0.2:1)*maxSpikeRate);
        if (row ==3)&&(col == 1)
            xlabel(ax, 'time (sec)');
            ylabel(ax, 'response');
        else
            set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        end
        axis(ax, 'square');
        text(ax,0.02, 180,sprintf('%4.2f c/deg', spatialFrequenciesCPD(sfIndex)), 'FontSize',16);
        
    end
    
    drawnow;
    exportFig = true;
    if (exportFig)
       plotlabOBJ.exportFig(hFig, 'pdf', sprintf('ResponseComponents_RGC_%2.0f_LMS_%0.2f_%0.2f_%0.2f', targetRGC, LMScontrast(1), LMScontrast(2), LMScontrast(3)), pwd());
    end
    setupPlotLab(-1);
end


function visualizeRetinalContrastSequence(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast)
        
    plotlabOBJ = setupPlotLab(0, 5, 6);
    hFig = figure(100); clf;
    
    midRow = round(size(theRetinalLMScontrastSequence,2)/2);
    colsNum = size(theRetinalLMScontrastSequence,3);
    margin = round(colsNum/4);
    colsToUse = margin:colsNum-margin;

    for frame = 1:1 %size(theRetinalLMScontrastSequence,1)
        
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,1)), 'Color', [1 0 0], 'LineWidth', 2); hold on;
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,2)), 'Color', [0 0.7 0], 'LineWidth', 2); 
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,3)), 'Color', [ 0 0 1], 'LineWidth', 2); 

        set(gca, 'YLim', [-0.11 0.11], 'XLim', [1 numel(colsToUse)], 'XTick', [], 'YTick', 0.1*[-1:0.5:1]);
        axis 'square'
        title(sprintf('%4.2f c/deg', gaborSpatialFrequencyCPD));
    end
    plotlab.offsetAxes(gca);
    drawnow;
    exportFig = true;
    if (exportFig)
       plotlabOBJ.exportFig(hFig, 'pdf', sprintf('RetinalContrastProfiles_%2.2fcpd_LMS_%0.2f_%0.2f_%0.2f', gaborSpatialFrequencyCPD, LMScontrast(1), LMScontrast(2), LMScontrast(3)), pwd());
    end
    setupPlotLab(-1);
end

function visualizeRetinalLMcontrastCorrelation(spatialFrequenciesCPD, roiLcontrast, roiMcontrast, LMScontrast)
    plotlabOBJ = setupPlotLab(0, 16, 12);
    hFig = figure(101); clf;
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 3, ...
        'colsNum', 4, ...
        'leftMargin', 0.06, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.03, ...
        'bottomMargin', 0.05, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    for sfIndex = 1:size(roiLcontrast,1)
        row = floor((sfIndex-1)/4)+1;
        col = mod(sfIndex-1,4)+1;
        ax = theAxesGrid{row,col};
        scatter(ax, roiLcontrast(sfIndex,1:4:end), roiMcontrast(sfIndex,1:4:end), '.');
        set(ax, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'XTick', 0.1*[-1:0.5:1], 'YTick', 0.1*[-1:0.5:1]);
        if (row ==3)&&(col == 1)
            xlabel(ax, 'retinal L-cone contrast');
            ylabel(ax, 'retinal M-cone contrast');
        else
            set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        end
        axis(ax, 'square');
        text(ax,-0.09, 0.085,sprintf('%4.2f c/deg', spatialFrequenciesCPD(sfIndex)), 'FontSize',16);
    end
    drawnow;
    exportFig = true;
    if (exportFig)
        plotlabOBJ.exportFig(hFig, 'pdf', sprintf('RetinalLMContrastCorrelationAtDifferentSFs_LMS_%0.2f_%0.2f_%0.2f', LMScontrast(1), LMScontrast(2), LMScontrast(3)), pwd());
    end
    setupPlotLab(-1);
    
end

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 

