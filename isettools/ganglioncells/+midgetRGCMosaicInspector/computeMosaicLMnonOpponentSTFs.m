function computeMosaicLMnonOpponentSTFs(mosaicCenterParams, mosaicSurroundParams, useParfor)
    
    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, mosaicSurroundParams.H1cellIndex);

    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');

    % Ask the user which optics position to use for the computation
    opticsPositionDegs = midgetRGCMosaicInspector.selectOpticsPosition(theMidgetRGCmosaic);

    % Generate the responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);
           
    viewingDistanceMeters = 4;
    stimulusPixelsNum = 512*2;
    coneContrasts = [1 1 0];
    deltaOri = 15;
    orientationsTested = 0:deltaOri:(180-deltaOri);
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 12 16 20 24 32 48 64];

    % Generate a presentation display with a desired resolution
    sceneFOVdegs = theMidgetRGCmosaic.inputConeMosaic.sizeDegs;
    retinalImageResolutionDegs = max(sceneFOVdegs)/stimulusPixelsNum;

    % At least 6 samples / period
    maxSF = 1/(2*3*retinalImageResolutionDegs);
    if (max(spatialFrequenciesTested) > maxSF)
        fprintf('Max SF examined (%2.2f c/deg) is too high for this FOV (%2.2f degs) and pixels num (%d). (SFmax: %2.2f c/deg)\n', ...
            max(spatialFrequenciesTested), max(sceneFOVdegs), stimulusPixelsNum, maxSF);
        idx = find(spatialFrequenciesTested <= maxSF);
        spatialFrequenciesTested = spatialFrequenciesTested(idx);
        if (maxSF > max(spatialFrequenciesTested))
            spatialFrequenciesTested(numel(spatialFrequenciesTested)+1) = maxSF;
        end

        fprintf('Will only measure the STF up to %2.2f c/deg.\n', max(spatialFrequenciesTested));
    end

    
    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            theMidgetRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', 0.75, ...
            'spatialFrequencyCPD', [], ...
            'orientationDegs', 0, ...
            'spatialPhaseIncrementDegs', 30, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', max(sceneFOVdegs), ...
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

    % Allocate memory
    stimParams.orientationDegs = 0;
    stimParams.spatialFrequencyCPD = spatialFrequenciesTested(1);
    [~, spatialPhasesDegs] = rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1);

    % Single precision responses
    theMidgetRGCMosaicResponses = zeros(...
        numel(orientationsTested), ...
        numel(spatialFrequenciesTested), ...
        numel(spatialPhasesDegs), ...
        rgcsNum, ...
        'single');
   
    disp('Allocated memory');
    
    % Go through all stimulus orientations
    for iOri = 1:numel(orientationsTested)
        stimParams.orientationDegs = orientationsTested(iOri);

        fprintf('Computing STF for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);

        for iFreq = 1:numel(spatialFrequenciesTested)

            theStimParams = stimParams;
            theStimParams.spatialFrequencyCPD = spatialFrequenciesTested(iFreq);
            
            % Generate spatial modulation patterns for each stimulus frame
            theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theStimParams);

            % Generate scenes for the different spatial phases
            [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, theStimParams, theDriftingGratingSpatialModulationPatterns, ...
                    'validateScenes', false);

            if (iOri == 1) && (iFreq == 1)
                % Do a compute just so we generate the optics
                theMidgetRGCmosaic.compute(...
                        theDriftingGratingFrameScenes{1}, ...
                        'nTrials', 1, ...
                        'theNullScene', theNullStimulusScene, ...
                        'withWavefronOpticsAtPositionDegs', opticsPositionDegs, ...
                        'normalizeConeResponsesWithRespectToNullScene', true);
            end


            % Allocate memory
            theFrameResponses = zeros(numel(spatialPhasesDegs), rgcsNum);

            if (useParfor)
            % Compute mRGCmosaic responses
                parfor iFrame = 1:numel(spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theScene = theDriftingGratingFrameScenes{iFrame};
    
                    % Compute the mosaic's response to this stimulus frame
                    r = theMidgetRGCmosaic.compute(...
                            theScene, ...
                            'nTrials', 1, ...
                            'theNullScene', theNullStimulusScene, ...
                            'withWavefronOpticsAtPositionDegs', opticsPositionDegs, ...
                            'normalizeConeResponsesWithRespectToNullScene', true);
    
                    theFrameResponses(iFrame,:) = r(1,1,:);
                end % iFrame
            else
                for iFrame = 1:numel(spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theScene = theDriftingGratingFrameScenes{iFrame};
    
                    % Compute the mosaic's response to this stimulus frame
                    r = theMidgetRGCmosaic.compute(...
                            theScene, ...
                            'nTrials', 1, ...
                            'theNullScene', theNullStimulusScene, ...
                            'withWavefronOpticsAtPositionDegs', opticsPositionDegs, ...
                            'normalizeConeResponsesWithRespectToNullScene', true);
    
                    theFrameResponses(iFrame,:) = r(1,1,:);
                end % iFrame
            end

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];

            theMidgetRGCMosaicResponses(iOri, iFreq,:,:) = single(theFrameResponses);

        end % iFreq
    end % iOri

    % Save all response data to disk
    save(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', 'opticsPositionDegs', '-v7.3');
end
