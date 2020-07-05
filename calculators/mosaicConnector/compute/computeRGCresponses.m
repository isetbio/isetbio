function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    presynapticSignal, spatialFrequenciesCPD, LMScontrast, saveDir)
    
    global LCONE_ID
    
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
           
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    LConeIndices = find(cmStruct.coneTypes == LCONE_ID);
    
    figure(55); clf;
    maxResponse = zeros(1, numel(spatialFrequenciesCPD));
    
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

        [centerResponses(sfIndex,:,:,:), surroundResponses(sfIndex,:,:,:)] = ...
            computeSubregionResponses(theConeMosaic, theMidgetRGCmosaic.centerWeights, theMidgetRGCmosaic.surroundWeights, thePresynapticResponses);

        fprintf('Done in %2.1f minutes\n', toc/60);
        
        % Display the mean presynaptic response
        showMeanMosaicResponseAsMovie = ~true;
        theMeanPresynapticResponses = squeeze(mean(thePresynapticResponses,1));
       
        if (showMeanMosaicResponseAsMovie)
            % Load the stimulus RGB sequence
            theStimulusRGBsequence = mFile.theStimulusRGBsequence;
            visualizeStimConeResponseMovie(responseTimeAxes, stimulusTimeAxis, theStimulusRGBsequence, theMeanPresynapticResponses, theConeMosaic);
        else
            
            subplot(3,5,sfIndex)
            imagesc(responseTimeAxis,1:numel(LConeIndices), theMeanPresynapticResponses(LConeIndices,:));
            set(gca, 'CLim', [0 600]);
            title('mean spatiotemporal response');
            colormap(gray);
            
            subplot(3,5,15)
            maxResponse(sfIndex) = max(max(theMeanPresynapticResponses(LConeIndices,:)));
            plot(spatialFrequenciesCPD(1:sfIndex), maxResponse(1:sfIndex), 'ks-');
            set(gca, 'XLim', [0.1 60], 'XScale', 'log');
           
        end
    end % sfIndex
    
    % Integrated center/surround responses
    integratedResponses = centerResponses - surroundResponses;
    
    % Average integrated responses over iterations
    integratedResponsesMean = mean(integratedResponses,2);
    integratedResponsesMean = integratedResponsesMean / max(abs(integratedResponsesMean(:)));
    
    centerResponsesMean = mean(centerResponses,2);
    surroundResponsesMean = mean(surroundResponses,2);
    m = max(centerResponsesMean(:));
    centerResponsesMean = centerResponsesMean/m;
    surroundResponsesMean = surroundResponsesMean/m;

    rgcsNum = size(integratedResponsesMean,3);
    figure(123); clf;
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        
        subplot(3,5,sfIndex);
        iRGC = 100;
        timeResponse = squeeze(integratedResponsesMean(sfIndex,1,iRGC,:));
        plot(1:length(timeResponse), timeResponse, 'gs-'); hold on;
        timeResponse = squeeze(centerResponsesMean(sfIndex,1,iRGC,:));
        plot(1:length(timeResponse), timeResponse, 'r-'); hold on;
        timeResponse = squeeze(surroundResponsesMean(sfIndex,1,iRGC,:));
        plot(1:length(timeResponse), timeResponse, 'b-'); hold on;
        
        set(gca, 'YLim', [-1 1]);

        
        title(sprintf('sf = %2.2f',spatialFrequenciesCPD(sfIndex)))
        for iRGC = 1:rgcsNum
            timeResponse = squeeze(integratedResponsesMean(sfIndex,1,iRGC,:));
            p = prctile(timeResponse,[15 85]);
            responseModulations(iRGC,sfIndex) = p(2)-p(1);
            prctiles(sfIndex,iRGC,:) = p;
        end
    end
    
    zLevels = [0.3 1];
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        visualizeRGCmosaicWithResponses(100+sfIndex, theConeMosaic, 'linear', ...
           responseTimeAxis, squeeze(integratedResponsesMean(sfIndex,1,:,:)), ...
           runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
           theMidgetRGCmosaic, zLevels, 'centers');    
    end
    
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, 'log', spatialFrequenciesCPD, responseModulations, ...
                runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                theMidgetRGCmosaic, zLevels, 'centers'); 
   
    
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
    
%     figure(222); clf;
%     subplot(1,2,1);
%     plot(t, centerIR, 'rs-'); hold on;
%     plot(t, surroundIR, 'bs-');
%     set(gca, 'XLim', [0 200/1000]);
    
    for instanceIndex = 1:instancesNum
        % All presynaptic spatiotemporal responses for this instance
        instancePresynapticResponse = squeeze(presynapticResponses(instanceIndex,:,:));
        for iRGC = 1:rgcsNum

            
            % The RGC's weights
            iRGCweightsC = (full(squeeze(weightsC(:,iRGC))))';
            iRGCweightsS = (full(squeeze(weightsS(:,iRGC))))';
            
%             if (iRGC == 100)
%                 figure(999)
%                 ax = subplot(2,2,1);
%                 theConeMosaic.renderActivationMap(ax, iRGCweightsC', 'signalRange', [0 max(iRGCweightsC)]);
%                 ax = subplot(2,2,2);
%                 theConeMosaic.renderActivationMap(ax, iRGCweightsS', 'signalRange', [0 max(iRGCweightsS)]);
%                 ax = subplot(2,2,3)
%                 theConeMosaic.renderActivationMap(ax, instancePresynapticResponse.*iRGCweightsC' );
%                  ax = subplot(2,2,4)
%                 theConeMosaic.renderActivationMap(ax, instancePresynapticResponse.*iRGCweightsS' );
%             end
            
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

%             [iRGC sum(iRGCweightsS)/sum(iRGCweightsC) mean(surroundR./centerR)]
             
%             subplot(1,2,2);
%             plot(t, centerR, 'r-');
%             hold on;
%             plot(t, surroundR, 'b-');
%             plot(t, centerR-surroundR, 'k-', 'LineWidth', 1.5);
%             pause
            
        end
    end
end