function preComputeConeMosaicLMnonOpponentSTFs(mosaicCenterParams, rfModelParams, opticsParams, useParfor)

    % Generate the frozen mosaic filename
    mosaicFileName = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);

    % Load the frozen midget RGC mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Ask the user which optics position to use for the computation
    opticsPositionDegs = midgetRGCMosaicInspector.selectOpticsPosition(theMidgetRGCmosaic);

    opticsParams
    pause
    % Generate the optics to use. Here we are using the wavefront measurements
    % (database, subject, eye) that were used to fit the midgetRGCmosaic RTVF objects,
    % but we can specify the retinal position of these optics, for example
    % some off-axis position, not at the center of the mosaic
    theMidgetRGCmosaic.generateOpticsAtPosition(opticsPositionDegs);
    theOI = theMidgetRGCmosaic.theCurrentOpticalImage;

    % Extract the input cone mosaic
    theInputConeMosaic = theMidgetRGCmosaic.inputConeMosaic;

    % Generate the input cone mosaics responses filename
    coneMosaicResponsesFileName = midgetRGCMosaicInspector.coneMosaicResponsesFileName(...
        mosaicFileName, opticsPositionDegs);

    % Generate components for running the STF mapping experiment
    [stimParams, thePresentationDisplay] = midgetRGCMosaicInspector.setupVisualSTFmappingExperiment(...
        theInputConeMosaic.sizeDegs, ...
        theInputConeMosaic.wave);

    % Allocate memory
    conesNum = size(theInputConeMosaic.coneRFpositionsDegs,1);

    % Single precision responses
    theConeMosaicResponses = zeros(...
        numel(stimParams.orientationsTested), ...
        numel(stimParams.spatialFrequenciesTested), ...
        numel(stimParams.spatialPhasesDegs), ...
        conesNum, ...
        'single');
   
    disp('Allocated memory');

    % Empty the noiseFreeAbsorptionsCountNull
    noiseFreeAbsorptionsCountNull = [];


    % Go through all stimulus orientations
    for iOri = 1:numel(stimParams.orientationsTested)
        stimParams.orientationDegs = stimParams.orientationsTested(iOri);

        fprintf('Computing cone mosaic STFs for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);

        for iFreq = 1:numel(stimParams.spatialFrequenciesTested)

            theCurrentStimParams = stimParams;
            theCurrentStimParams.spatialFrequencyCPD = stimParams.spatialFrequenciesTested(iFreq);
            
            % Generate spatial modulation patterns for each stimulus frame
            theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theCurrentStimParams);

            % Generate scenes for the different spatial phases
            [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                    thePresentationDisplay, theCurrentStimParams, theDriftingGratingSpatialModulationPatterns, ...
                    'validateScenes', false);


            % Compute the cone mosaic responses to the null scene OI
            if (isempty(noiseFreeAbsorptionsCountNull))
                % Compute the optical image of the null scene
                theOI  = oiCompute(theNullStimulusScene, theOI);

                % Compute the cone mosaic responses
                coneMosaicNullResponses = theInputConeMosaic.compute(theOI, 'nTrials', 1);
            end



            % Allocate memory
            theFrameResponses = zeros(numel(theCurrentStimParams.spatialPhasesDegs), conesNum, 'single');

            % Compute the input cone mosaic responses
            if (useParfor)
                parfor iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    theCurrentOI = oiCompute(theFrameScene, theOI);

                    % Compute the cone mosaic responses
                    [noiseFreeAbsorptionsCount, noisyAbsorptionsCountInstances, noiseFreePhotoCurrents, noisyPhotocurrentInstances, ~] = ...
                        theInputConeMosaic.compute(theCurrentOI, 'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));
                end
            else
                for iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    theCurrentOI  = oiCompute(theFrameScene, theOI);

                    % Compute the cone mosaic responses
                    [noiseFreeAbsorptionsCount, noisyAbsorptionsCountInstances, photoCurrents, ~, ~] = ...
                        theInputConeMosaic.compute(theCurrentOI, 'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));
                end
            end

            % Save theFrameResponses
            theConeMosaicResponses(iOri, iFreq,:,:) = theFrameResponses;

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
    save(coneMosaicResponsesFileName, 'theConeMosaicResponses', 'coneMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', 'opticsPositionDegs', '-v7.3');

    fprintf('Saved computed cone mosaic responses to %s\n', coneMosaicResponsesFileName);
end

