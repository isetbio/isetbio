%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - how to select one of the pre-computed midget RGC mosaics,
%        - how to select noise source (cone-only or intrinscic mRGC noise only)
%        - how to compute responses to a static (single-frame) stimulus
%        - how to visualize different aspects of the mRGCMosaic
%        - how to visualize its response
%


% History:
%    10/25/23  NPC  Wrote it.

function t_mRGCMosaicBasic
    %% Close all figures
    close all;

    %% Display available mRGCMosaics
    rgcMosaicType = 'ONcenterMidgetRGC';
    mRGCMosaic.availableComputeReadyMosaics(rgcMosaicType);

    %% Specify the desired eccentricity of the precomputed mRGC mosaic
    % Choose the x-eccentricity from one of the available mosaics displayed above
    % (e.g., -16.0 to load the mosaic 'mRGCMosaicEcDegs(-10.0_0.0)_SizeDegs(6.0_3.0)...'
    horizontalEccDegs = input('Enter mRGCMosaic''s horizontal eccentricity: ');

    %% Load precomputed mRGCMosaic
    theMRGCMosaic = MosaicPoolingOptimizer.loadPreComputedMRGCMosaic(horizontalEccDegs);

    %% Set the presentation display on which the stimulus must be realized
    presentationDisplay = generateCRTDisplay();

    %% Set stimulus parameters
    sinusoidalStimulusParams = struct(...
        'sizeDegs', 2.0, ...
        'positionDegs',[horizontalEccDegs 0.], ...
        'spatialFrequencyCyclesPerDeg', 10.0, ...
        'orientationDegs', 0, ...               
        'phaseDegs', 0, ...  
        'contrast', 0.9, ...
        'meanLuminanceCdPerM2', 40, ...
        'isWindowed', false, ...
        'pixelsAlongWidthDim', 1024, ...       
        'pixelsAlongHeightDim', 1024);

    %% Generate a sinusoidal stimulus scene and the background scene
    [theStimulusScene, theBackgroundScene] = generateSinusoidalStimulusAndBackgroundScene(...
        presentationDisplay, sinusoidalStimulusParams);

    %% Retrieve the native optics
    % These are the optics under which the connections from the input cone
    % mosaic to the mRGC mosaic were optimized
    theOI = theMRGCMosaic.theNativeOptics;

    %% Compute the retinal image of the stimulus and background scenes
    theStimulusRetinalImage = oiCompute(theStimulusScene, theOI);
    theBackgroundRetinalImage = oiCompute(theBackgroundScene, theOI);

    %% Set the integration time of the input cone mosaic to 100 msec
    theMRGCMosaic.inputConeMosaic.integrationTime = 100/1000;

    %% Set the input cone mosaic noise flag to none
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'none';
    
    %% Compute the noise-free background activation of the input cone mosaic
    theConeMosaicBackgroundNoiseFreeActivation = theMRGCMosaic.inputConeMosaic.compute(...
            theBackgroundRetinalImage, ...
            'opticalImagePositionDegs', sinusoidalStimulusParams.positionDegs);

    %% Set the input cone mosaic noise flag to random to generate noisy response instances
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'random';

    %% Number of noisy response instances to compute
    noisyInstancesNum = 128;

    %% Compute the noise-free response and noisyInstancesNum of noisy response instances of the input cone mosaic
    [theConeMosaicNoiseFreeResponse,  ...
     theConeMosaicNoisyResponseInstances, ~, ~, ...
     theConeMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.inputConeMosaic.compute(...
            theStimulusRetinalImage, ...
            'opticalImagePositionDegs', sinusoidalStimulusParams.positionDegs, ...
            'nTrials', noisyInstancesNum);

    %% Visualize the spatial relationship between stimulus and MRGCMmosaic
    % This function must be called after the cMosaic.compute() method in
    % order to display the OI at the correct location within the mosaic
    coVisualizeRetinalStimulusConeAndMRGCmosaic(1, theMRGCMosaic, theStimulusRetinalImage);

    %% Define a function handle to convert excitations to modulations
    excitationsToModulations = @(e) (bsxfun(@times, (bsxfun(@minus, e, theConeMosaicBackgroundNoiseFreeActivation)), 1./theConeMosaicBackgroundNoiseFreeActivation));

    %% Set the mRGCMosaic noise flag to none
    theMRGCMosaic.noiseFlag = 'none';

    %% Compute the noise-free response of the mRGC mosaic. No cone noise, no mRGC noise
    noiseFreeResponse = theMRGCMosaic.compute( ...
              excitationsToModulations(theConeMosaicNoiseFreeResponse), ...
              theConeMosaicResponseTemporalSupportSeconds);

    %% Set the mRGCMosaic noise flag to random
    theMRGCMosaic.noiseFlag = 'random';
 
    %% Set the intrinsic noise of the mRGCMosaic to 0
    theMRGCMosaic.vMembraneGaussianNoiseSigma = 0.0;

    %% Compute noisyInstancesNum of mRGC mosaic responses with the cone mosaic noise being the only noise source
    % To do so, we compute using the input cone mosaic noisy modulation response instances
    [~, noisyResponseInstancesConeNoiseOnly] = theMRGCMosaic.compute( ...
             excitationsToModulations(theConeMosaicNoisyResponseInstances), ...
             theConeMosaicResponseTemporalSupportSeconds, ...
             'nTrials', noisyInstancesNum);

    %% Set the intrinsic noise of the mRGCMosaic to the examined level
    theMRGCMosaic.vMembraneGaussianNoiseSigma = 0.05;

    %% Compute noisyInstancesNum of mRGC mosaic responses with the mRGCMosaic intrinsic noise being the only noise source
    % To do so, we compute using the input cone mosaic noise-free modulation response 
    [~, noisyResponseInstancesIntrinsicMRGCnoiseOnly] = theMRGCMosaic.compute( ...
             excitationsToModulations(theConeMosaicNoiseFreeResponse), ...
             theConeMosaicResponseTemporalSupportSeconds, ...
             'nTrials', noisyInstancesNum);

    %% Find indices of RGCs along the y = y-stimulus position
    targetYdegs = sinusoidalStimulusParams.positionDegs(2);
    minConeSeparation = 0;
    minRGCSeparation = 0;
    [~, mRGCIndices] = coneAndMRGCindicesAlongDesiredYposition(...
        theMRGCMosaic, targetYdegs, ...
        minConeSeparation, minRGCSeparation);

    %% Visualize response components
    generatePDFs = ~true;
    visualizeMRGMosaicResponseComponents(theMRGCMosaic, mRGCIndices, ...
        noiseFreeResponse, ...
        noisyResponseInstancesConeNoiseOnly, ...
        noisyResponseInstancesIntrinsicMRGCnoiseOnly, generatePDFs);
end