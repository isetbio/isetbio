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
           
    % Generate components for running the STF mapping experiment
    [stimParams, thePresentationDisplay] = midgetRGCMosaicInspector.setupVisualSTFmappingExperiment(...
        theMidgetRGCmosaic.inputConeMosaic.sizeDegs, ...
        theMidgetRGCmosaic.inputConeMosaic.wave);

    % Allocate memory
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1);

    % Single precision responses
    theMidgetRGCMosaicResponses = zeros(...
        numel(stimParams.orientationsTested), ...
        numel(stimParams.spatialFrequenciesTested), ...
        numel(stimParams.spatialPhasesDegs), ...
        rgcsNum, ...
        'single');
   
    disp('Allocated memory');
    
    % Initialize theInputDataStruct for computation with input scenes
    theInputDataStruct = midgetRGCMosaic.inputDataStruct(midgetRGCMosaic.SCENE_COMPUTE_INPUT_DATA_TYPE);


    % Go through all stimulus orientations
    for iOri = 1:numel(stimParams.orientationsTested)
        stimParams.orientationDegs = stimParams.orientationsTested(iOri);

        fprintf('Computing midget RGC mosaic STFs for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);

        for iFreq = 1:numel(stimParams.spatialFrequenciesTested)

            theCurrentStimParams = stimParams;
            theCurrentStimParams.spatialFrequencyCPD = stimParams.spatialFrequenciesTested(iFreq);
            
            % Generate spatial modulation patterns for each stimulus frame
            theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theCurrentStimParams);

            % Generate scenes for the different spatial phases
            [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    thePresentationDisplay, theCurrentStimParams, theDriftingGratingSpatialModulationPatterns, ...
                    'validateScenes', false);


            if (iOri == 1) && (iFreq == 1)
                theCurrentInputDataStruct = theInputDataStruct;
                theCurrentInputDataStruct.theTestScene = theDriftingGratingFrameScenes{1}; 
                theCurrentInputDataStruct.theNullScene = theNullStimulusScene;
                theCurrentInputDataStruct.wavefrontOpticsPositionDegs = theMidgetRGCmosaic.eccentricityDegs;
                theCurrentInputDataStruct.opticalImagePositionDegs = 'mosaic-centered';
                theCurrentInputDataStruct.normalizeConeResponsesWithRespectToNullScene = true;
            
                % Do a compute just so we generate the optics
                theMidgetRGCmosaic.compute(...
                    theCurrentInputDataStruct, ...
                    'nTrials', 1);
            end


            % Allocate memory
            theFrameResponses = zeros(numel(theCurrentStimParams.spatialPhasesDegs), rgcsNum);

            % Compute mRGCmosaic responses
            if (useParfor)

                % Compute each frame separately
                parfor iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};
    
                    % Update theInputDataStruct with theFrameScene
                    theCurrentInputDataStruct = theInputDataStruct;
                    theCurrentInputDataStruct.theTestScene = theFrameScene; 
                    theCurrentInputDataStruct.theNullScene = theNullStimulusScene;
                    theCurrentInputDataStruct.wavefrontOpticsPositionDegs = theMidgetRGCmosaic.eccentricityDegs;
                    theCurrentInputDataStruct.opticalImagePositionDegs = 'mosaic-centered';
                    theCurrentInputDataStruct.normalizeConeResponsesWithRespectToNullScene = true;

                    % Compute the mosaic's response to this stimulus frame
                    r =  theMidgetRGCmosaic.compute(...
                            theCurrentInputDataStruct, ...
                            'nTrials', 1);
    
                    % Insert it in theFrameResponses
                    theFrameResponses(iFrame,:) = r(1,1,:);
                end % iFrame
            else
                for iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};
    
                    % Update theInputDataStruct with theFrameScene
                    theCurrentInputDataStruct = theInputDataStruct;
                    theCurrentInputDataStruct.theTestScene = theFrameScene; 
                    theCurrentInputDataStruct.theNullScene = theNullStimulusScene;
                    theCurrentInputDataStruct.wavefrontOpticsPositionDegs = theMidgetRGCMosaic.eccentricityDegs;
                    theCurrentInputDataStruct.opticalImagePositionDegs = 'mosaic-centered';
                    theCurrentInputDataStruct.normalizeConeResponsesWithRespectToNullScene = true;

                    % Compute the mosaic's response to this stimulus frame
                    r =  theMidgetRGCmosaic.compute(...
                            theCurrentInputDataStruct, ...
                            'nTrials', 1);
    
                    theFrameResponses(iFrame,:) = r(1,1,:);
                end % iFrame
            end

            % Save theFrameResponses
            theMidgetRGCMosaicResponses(iOri, iFreq,:,:) = single(theFrameResponses);

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];
        end % iFreq
    end % iOri

    % Save response data to disk
    orientationsTested = stimParams.orientationsTested;
    spatialFrequenciesTested = stimParams.spatialFrequenciesTested;
    coneContrasts = stimParams.coneContrasts;
    spatialPhasesDegs = stimParams.spatialPhasesDegs;
    save(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', 'opticsPositionDegs', '-v7.3');
end
