function computeMosaicLMnonOpponentSTFs(mosaicCenterParams, mosaicSurroundParams)
    
    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.generateMosaicFileName(...
        mosaicCenterParams);
   
    % Generate responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileNameForMosaicFileName(...
        mosaicFileName, mosaicSurroundParams.H1cellIndex);
           
    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

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
    theMidgetRGCMosaicResponses = ...
        zeros(numel(orientationsTested), numel(spatialFrequenciesTested), numel(spatialPhasesDegs), rgcsNum);
   
    disp('Allocated memory');
    
    % No optics. We will get it back during the first call to
    % mRGCMosaic.compute()
    theOptics = [];

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

            % Allocate memory
            theFrameResponses = zeros(numel(spatialPhasesDegs), rgcsNum);

            % Compute mRGCmosaic responses
            for iFrame = 1:numel(spatialPhasesDegs)

                fprintf('Computing mRGC mosaic response to frame (%d/%d) of the %2.2f c/deg stimulus.\n', ...
                    iFrame, numel(spatialPhasesDegs), theStimParams.spatialFrequencyCPD);

                % Get scene corresponding to this stimulus frame
                theScene = theDriftingGratingFrameScenes{iFrame};

                if (isempty(theOptics))
                    % Compute the mosaic's response to this stimulus frame
                    % and also retrieve the optics so we can pass it along
                    % in subsequent calls (avoid recomputing it)
                    [r,~,theConeMosaicActivation, theOptics] = theMidgetRGCmosaic.compute(...
                        theScene, ...
                        'nTrials', 1, ...
                        'theNullScene', theNullStimulusScene, ...
                        'normalizeConeResponsesWithRespectToNullScene', true);
                else
                    % Compute the mosaic's response to this stimulus frame
                    % using the returned optics
                    [r,~,theConeMosaicActivation] = theMidgetRGCmosaic.compute(...
                        theScene, ...
                        'nTrials', 1, ...
                        'theNullScene', theNullStimulusScene, ...
                        'withOptics', theOptics, ...
                        'normalizeConeResponsesWithRespectToNullScene', true);
                end


                % Store the mosaic's responses
                theFrameResponses(iFrame,:) = r;

            end % iFrame

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];

            theMidgetRGCMosaicResponses(iOri, iFreq,:,:) = theFrameResponses;

        end % iFreq
    end % iOri

    % Save all response data to disk
    save(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', '-v7.3');
end
