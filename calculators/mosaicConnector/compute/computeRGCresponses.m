function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    presynapticSignal, spatialFrequenciesCPD, LMScontrast, saveDir)
    
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
        
        % Compute mean over all instances integrated response
        integratedResponsesMean(sfIndex,:,:) = squeeze(mean(cR-sR, 1));
        
        % Compute modulation of the response at the fundamental temporal frequency
        fundamentalTuning(:, sfIndex) = computeFundamentalResponseModulation(squeeze(integratedResponsesMean(sfIndex,:,:)));
        
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
    
    
    iRGC = 100;
    figNo = 123;
    visualizeResponseComponentsForSingleRGC(figNo, iRGC, responseTimeAxis, centerResponseInstances, surroundResponseInstances, spatialFrequenciesCPD);  
    
    % Visualize the temporal response of each RGC at the RGC's location
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        zLevels = [0.3 1];
        visualizeRGCmosaicWithResponses(100+sfIndex, theConeMosaic, 'linear', ...
           responseTimeAxis, squeeze(integratedResponsesMean(sfIndex,:,:)), ...
           runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
           theMidgetRGCmosaic, zLevels, 'centers');    
    end
    
    % Visualize the spatial frequency tuning of each RGC at the RGC's location
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, 'log', spatialFrequenciesCPD, fundamentalTuning, ...
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

function fundamentalTuning = computeFundamentalResponseModulation(responses)
        % Compute response modulations
        rgcsNum = size(responses,1);
        fundamentalTuning = zeros(rgcsNum,1);
        for iRGC = 1:rgcsNum
            timeResponse = squeeze(responses(iRGC,:));
            p = prctile(timeResponse,[15 85]);
            fundamentalTuning(iRGC) = p(2)-p(1);
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
        
function  visualizeResponseComponentsForSingleRGC(figNo, iRGC, responseTimeAxis, centerResponses, surroundResponses, spatialFrequenciesCPD)
    centerResponsesMean = squeeze(mean(centerResponses,2));
    surroundResponsesMean = squeeze(mean(surroundResponses,2));
    m = max(abs(centerResponsesMean(:)));

    figure(figNo); clf;
    for sfIndex = 1:numel(spatialFrequenciesCPD)    
        subplot(3,5,sfIndex);
        centerResponses = squeeze(centerResponsesMean(sfIndex,iRGC,:));
        plot(responseTimeAxis, centerResponses, 'r-'); hold on;
        surroundResponses = squeeze(surroundResponsesMean(sfIndex,iRGC,:));
        plot(responseTimeAxis, surroundResponses, 'b-'); hold on;
        set(gca, 'YLim', [-m m]);
    end
end
