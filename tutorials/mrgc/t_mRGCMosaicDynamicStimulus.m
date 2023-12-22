%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - how to select one of the pre-computed midget RGC mosaics,
%        - how to compute responses to a dynamic stimulus, and
%        - how to visualize its dynamic response
%


% History:
%    10/25/23  NPC  Wrote it.

function t_mRGCMosaicDynamicStimulus
    %% Close all figures
    close all;

    % Configure a conservative parpool manager. This gives at least 8 GB RAM/core
    ASPPManager = AppleSiliconParPoolManager('conservative');

    % Control saving of figures.  We don't want tutorials
    % saving things into the isetbio source tree.
    saveFigures = true;
    figureDir = figExporter.figureDir(mfilename, saveFigures);
    
    %% Display available mRGCMosaics
    rgcMosaicType = 'ONcenterMidgetRGC';
    mRGCMosaic.availableComputeReadyMosaics(rgcMosaicType);

    %% Specify the desired eccentricity of the precomputed mRGC mosaic
    % Choose the x-eccentricity from one of the available mosaics,
    % displayed above
    % (e.g., -16.0 to load the mosaic 'mRGCMosaicEcDegs(-10.0_0.0)_SizeDegs(6.0_3.0)...'
    horizontalEccDegs = 2.5; %input('Enter mRGCMosaic''s horizontal eccentricity: ');

    %% Load the precomputed mRGCMosaic
    theMRGCMosaic = MosaicPoolingOptimizer.loadPreComputedMRGCMosaic(horizontalEccDegs);

    %% Examined level of the mRGCMosaic intrinsic noise
    intrinsicMRGCnoiseSigma = 0.20;

    %% Presentation display params
    viewingDistanceMeters = 1.0;
    wavelengthSupport = theMRGCMosaic.inputConeMosaic.wave;

    % Stimulus resolution so that cone aperture blur will have an observable effect
    pixelSizeDegs = min(theMRGCMosaic.inputConeMosaic.coneApertureDiametersDegs)/5;
      
    %% Generate the presentation display on which the stimulus must be realized
    presentationDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            wavelengthSupport, pixelSizeDegs, ...
            viewingDistanceMeters);

    %% Set drifting sinusoidal stimulus parameters
    %
    % These are set with a short duration so this runs reasonably
    % fast.  Lengthen durationSecs (e.g to 1) if you want to
    % see more of the drift.  But beware, it will take while
    % and might fill your computer's memory and crash if you do
    % so.
    driftingSinusoidalStimulusParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', [1 1 1], ...
            'contrast', 0.5, ...
            'orientationDegs', 90, ...
            'spatialFrequencyCPD', 2.0, ...
            'spatialPhaseIncrementDegs', 45, ...
            'temporalFrequencyHz', 2.0, ...
            'durationSeconds', 0.25, ...
            'temporalEnvelopeTau', 0.5/3, ...
            'pixelSizeDegs', pixelSizeDegs, ...
            'stimSizeDegs', 4.0, ...
            'positionDegs', theMRGCMosaic.eccentricityDegs, ...
            'wavelengthSupport', displayGet(presentationDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(presentationDisplay, 'viewing distance') ...
            );

    %% Generate spatial modulation patterns for each stimulus frame
    [driftingGratingSpatialModulationPatterns, ...
     driftingSinusoidalStimulusParams.spatialPhasesDegs, ...
     theStimulusTemporalSupportSeconds, theStimulusTemporalRamp] = ...
        rfMappingStimulusGenerator.driftingGratingFrames(driftingSinusoidalStimulusParams);


    %% Generate a sequence of scenes representing a drifting sinusoidal stimulus and the background scene
    [theDriftingGratingSceneSequence, theBackgroundScene] = ...
        rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
           presentationDisplay, driftingSinusoidalStimulusParams, driftingGratingSpatialModulationPatterns, ...
           'validateScenes', ~true);

    %% Retrieve the native optics
    % These are the optics under which we optimized connections from the input cone
    % mosaic to the mRGC mosaic 
    theOI = theMRGCMosaic.theNativeOptics;

    %% Compute the retinal image of the background scene
    theBackgroundRetinalImage = oiCompute(theOI, theBackgroundScene);

    %% Compute the retinal image sequence corresponding to the different frames of the drifting grating
    framesNum = numel(theDriftingGratingSceneSequence);
    theListOfRetinalImages = cell(1, framesNum);
    for iFrame = 1:framesNum
        fprintf('Generating retinal image sequence for frame %d of %d\n', iFrame, framesNum);
        theListOfRetinalImages{iFrame} = ...
            oiCompute(theOI,theDriftingGratingSceneSequence{iFrame});
    end

    %% Generate an @oiSequence object from the list of computed optical images
    theStimulusRetinalOISequence = oiArbitrarySequence(...
        theListOfRetinalImages, theStimulusTemporalSupportSeconds);

    %% Save some RAM
    clear 'theListOfRetinalImages'

    %% Set the input cone mosaic integration time to the duration of a single stimulus frame
    theMRGCMosaic.inputConeMosaic.integrationTime = ...
        theStimulusRetinalOISequence.timeAxis(2) - theStimulusRetinalOISequence.timeAxis(1);

    %% Set the input cone mosaic noise flag to none
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'none';

    %% Compute the spatiotemporal noise-free cone-mosaic activation
    [theConeMosaicNoiseFreeSpatiotemporalResponse, ~, ~, ~, ...
     theConeMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.inputConeMosaic.compute(...
            theStimulusRetinalOISequence, ...
            'opticalImagePositionDegs', driftingSinusoidalStimulusParams.positionDegs);

   
    %% Compute the noise-free background activation of the input cone mosaic
    theConeMosaicNoiseFreeBackgroundResponse = theMRGCMosaic.inputConeMosaic.compute(...
            theBackgroundRetinalImage, ...
            'opticalImagePositionDegs', driftingSinusoidalStimulusParams.positionDegs);

    %% Visualize the spatial relationship between stimulus and MRGCMmosaic
    % This function must be called after the cMosaic.compute() method in
    % order to display the OI at the correct location within the mosaic
    [~,maxModulationFrame] = max(theStimulusTemporalRamp(:));
    figNo = 1;
    hFig = coVisualizeRetinalStimulusConeAndMRGCmosaic(figNo, theMRGCMosaic, theStimulusRetinalOISequence.frameAtIndex(maxModulationFrame));

    % Save figure
    pdfFileName = fullfile(figureDir, 'mRGCmosaic_ConeMosaic_OpticalImage_Combo.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    %% Save some RAM
    clear 'theStimulusRetinalOISequence'

    %% Define a function handle to convert excitations to modulations
    excitationsToModulations = @(e) (bsxfun(@times, (bsxfun(@minus, e, theConeMosaicNoiseFreeBackgroundResponse)), 1./theConeMosaicNoiseFreeBackgroundResponse));

    %% Number of noisy response instances to compute
    noisyInstancesNum = 4;

    %% Set the mRGCMosaic noise flag
    theMRGCMosaic.noiseFlag = 'random';

    %% Set the intrinsic noise of the mRGCMosaic to the examined level
    theMRGCMosaic.vMembraneGaussianNoiseSigma = intrinsicMRGCnoiseSigma;

    %% Compute the mRGC spatiotemporal response
    [mRGCMosaicNoiseFreeSpatiotemporalResponse, ...
     mRGCMosaicNoisySpatiotemporalResponseInstances, ...
     mRGCMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
             excitationsToModulations(theConeMosaicNoiseFreeSpatiotemporalResponse), ...
             theConeMosaicResponseTemporalSupportSeconds, ...
             'nTrials', noisyInstancesNum);


    %% Find indices of RGCs along the y-stimulus position
    targetYdegs = driftingSinusoidalStimulusParams.positionDegs(2);
    minConeSeparation = 0;
    minRGCSeparation = 0.05;
    [~, mRGCIndices] = coneAndMRGCindicesAlongDesiredYposition(...
        theMRGCMosaic, targetYdegs, ...
        minConeSeparation, minRGCSeparation);

    % Export the video
    exportVideo = false;
    videofFileName = fullfile(isetbioRootPath, 'local', 'mRGCmosaicDynamicActivationVideo');

    figNo =  figNo + 1;
    visualizeDynamicMRGCMosaicResponse(figNo, theMRGCMosaic, mRGCIndices, ...
        mRGCMosaicNoiseFreeSpatiotemporalResponse, ...
        mRGCMosaicNoisySpatiotemporalResponseInstances,  ...
        mRGCMosaicResponseTemporalSupportSeconds, ...
        exportVideo, videofFileName);

    % Restore previous parpool config
    ASPPManager.restoreLastParpoolSize();
end