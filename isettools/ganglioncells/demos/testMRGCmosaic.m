function testMRGCmosaic

    recomputeReponses = true;
    if (recomputeReponses)
        
        % Instantiate a [0.3 x 0.3] deg wide midget RGC mosaic positioned 
        % at an eccentricity of [-1 0] degs , in the right eye
        mosaicEccDegs = [7 0];
        mosaicSizeDegs = [1 1];
        whichEye = 'right';

        % Specify spatial and temporal parameters of the drifting grating sequence
        presentationMode = 'drifted';

        % Temporal frequency in Hz
        temporalFrequencyHz = 10;    

        % Stimulus duration in seconds
        stimulusDurationSeconds = 0.5;

        % How motion is sampled, 45 degs = 8 spatial phases/period
        spatialPhaseAdvanceDegs = 45;    

        frameDurationSeconds = 1.0/(360/spatialPhaseAdvanceDegs*temporalFrequencyHz)
        pause
        
        % Spatial frequency 
        spatialFrequencyCPD = 1.0;

        % Specify chromatic direction
        chromaDir = [1.0, 1.0, 1.0]';

        % Instantiate a scene engine for drifting gratings
        driftingGratingSceneEngine = createGratingScene(chromaDir, spatialFrequencyCPD, ...
            'duration', stimulusDurationSeconds, ...
            'temporalFrequencyHz', temporalFrequencyHz, ...
            'spatialPhaseAdvanceDegs', spatialPhaseAdvanceDegs, ...
            'fovDegs', max(mosaicSizeDegs)*2, ...
            'spatialEnvelope', 'square', ...
            'presentationMode', presentationMode ...
            );

        % Generate a drifting grating sequence with 50% contrast
        sceneContrast = 0.5;
        [theDriftingGratingSequence, theStimulusTemporalSupportSeconds] = driftingGratingSceneEngine.compute(sceneContrast);

        % Visualize the drifting sequence
        driftingGratingSceneEngine.visualizeSceneSequence(theDriftingGratingSequence, theStimulusTemporalSupportSeconds);

        % Generate default optics
        theOptics = oiCreate();

        % Compute the sequence of optical images corresponding to the sequence of scenes
        framesNum = numel(theDriftingGratingSequence);
        theListOfOpticalImages = cell(1, framesNum);
        for frame = 1:framesNum
            theListOfOpticalImages{frame} = oiCompute(theDriftingGratingSequence{frame}, theOptics);
        end

        % Generate an @oiSequence object containing the list of computed optical images
        theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theStimulusTemporalSupportSeconds);

        % Generate a midget RGC mosaic
        theMidgetRGCmosaic = mRGCmosaic(mosaicEccDegs, mosaicSizeDegs, whichEye, ...
            'coneSpecificityLevel', 100, ...
            'viewTesselationMaps', ~true, ...
            'coneMosaicIntegrationTime', frameDurationSeconds);

        % Visualize mRGC positions together with imported ecc-varying cone positions
        % and together with cone positions in the equivalent employed reg-hex mosaic
        theMidgetRGCmosaic.visualizeInputPositions();

        % Visualize connections of cones to mRGC RF centers
        theMidgetRGCmosaic.visualizeConeMosaicTesselation('degrees');

        %theMidgetRGCmosaic.visualizeSynthesizedParams();
        %theMidgetRGCmosaic.visualizeConeWeights();


        % Compute responses
        instancesNum = 2;

        % Retrieve the coneMosaic that provides input to the mRGC mosaic
        theConeMosaic = theMidgetRGCmosaic.inputConeMosaic;

        % Zero eye movements
        eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
        emPaths = zeros(instancesNum, eyeMovementsNum, 2);

        % Compute the cone mosaic responses
        theConeMosaic.noiseFlag = 'none';
        theMidgetRGCmosaic.noiseFlag = 'none';
        
        theConeMosaicResponses = theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...   % the emPaths
            'currentFlag', false ...  % no photocurrent
        );

        % Compute the mRGC mosaic responses
        [theMRGCMosaicResponses, temporalSupportSeconds] = theMidgetRGCmosaic.compute( ...
                    theConeMosaicResponses, theConeMosaic.timeAxis);


        save('tmp.mat', ...
            'theMidgetRGCmosaic', ...
            'temporalSupportSeconds', ...
            'theMRGCMosaicResponses', ...
            'theStimulusTemporalSupportSeconds', ...
            'theDriftingGratingSequence', ...
            'theOIsequence', ...
            'theConeMosaicResponses', ...
            '-v7.3');
    else
        load('tmp.mat', ...
            'theMidgetRGCmosaic', ...
            'temporalSupportSeconds', ...
            'theMRGCMosaicResponses', ...
            'theStimulusTemporalSupportSeconds', ...
            'theDriftingGratingSequence', ...
            'theOIsequence', ...
            'theConeMosaicResponses');
    end
    
    
    % Visualize responses
    theMidgetRGCmosaic.visualizeResponses(temporalSupportSeconds, ...
        theMRGCMosaicResponses, ...
        'stimulusTemporalSupportSeconds', theStimulusTemporalSupportSeconds,...
        'stimulusSceneSequence', theDriftingGratingSequence, ...
        'opticalSequence', theOIsequence, ...
        'coneMosaicResponse', theConeMosaicResponses);
    
   
    
end

