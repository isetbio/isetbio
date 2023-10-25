%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - how to select one of the pre-computed midget RGC mosaics,
%        - how to compute responses to a static (single-frame) stimulus, and
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
    % Choose the x-eccentricity from one of the available mosaics,
    % displayed above
    % (e.g., -16.0 to load the mosaic 'mRGCMosaicEcDegs(-10.0_0.0)_SizeDegs(6.0_3.0)...'
    horizontalEccDegs = input('Enter mRGCMosaic''s horizontal eccentriciy: ');

    %% Load the precomputed mRGCMosaic
    theMRGCMosaic = MosaicPoolingOptimizer.loadPreComputedMRGCMosaic(horizontalEccDegs);

    %% Examined level of the mRGCMosaic intrinsic noise
    intrinsicMRGCnoiseSigma = 0.20;

    %% Set the presentation display on which the stimulus must be realized
    presentationDisplay = generateCRTDisplay();

    %% Set stimulus parameters
    stimulusSizeDegs = 10.0;        
    stimulusPositionDegs = [7 0];
    stimulusSpatialFrequencyCPD = 2.0;
    stimulusIsWindowed = false;
    stimulusContrast = 0.9;

    %% Generate a sinusoidal stimulus scene and the background scene
    [theStimulusScene, theBackgroundScene] = generateSinusoidalStimulusAndBackgroundScene(...
        presentationDisplay, stimulusSizeDegs, stimulusSpatialFrequencyCPD, ...
        stimulusContrast, stimulusIsWindowed);

    %% Retrieve the native optics
    % These are the optics under which we optimized connections from the input cone
    % mosaic to the mRGC mosaic 
    theOI = theMRGCMosaic.theNativeOptics;

    %% Compute the retinal image of the stimulus and background scenes
    theStimulusRetinalImage = oiCompute(theStimulusScene, theOI);
    theBackgroundRetinalImage = oiCompute(theBackgroundScene, theOI);

    %% Visualize the spatial relationship between stimulus and MRGCMmosaic
    visualizeSpatialRelationshipBetweenStimulusAndMRGCmosaic(1, theMRGCMosaic, theStimulusRetinalImage);

    %% Number of noisy response instances to compute
    noisyInstancesNum = 128;

    %% Set the input cone mosaic noise flag
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'random';

    %% Compute the noise-free response and noisyInstancesNum of noisy response instances of the input cone mosaic
    [theConeMosaicNoiseFreeResponse,  ...
     theConeMosaicNoisyResponseInstances, ~, ~, ...
     theConeMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.inputConeMosaic.compute(...
            theStimulusRetinalImage, ...
            'opticalImagePositionDegs', stimulusPositionDegs, ...
            'nTrials', noisyInstancesNum);

    %% Compute the background activation of the input cone mosaic
    theConeMosaicBackgroundActivation = theMRGCMosaic.inputConeMosaic.compute(...
            theBackgroundRetinalImage, ...
            'opticalImagePositionDegs', stimulusPositionDegs);

    %% Define a Function handle to convert excitations to modulations
    excitationsToModulations = @(e) (bsxfun(@times, (bsxfun(@minus, e, theConeMosaicBackgroundActivation)), 1./theConeMosaicBackgroundActivation));

    %% Compute the noise-free response of the mRGC mosaic. No cone noise, no mRGC noise
    theMRGCMosaic.noiseFlag = 'none';
    noiseFreeResponse = theMRGCMosaic.compute( ...
              excitationsToModulations(theConeMosaicNoiseFreeResponse), ...
              theConeMosaicResponseTemporalSupportSeconds);

    %% Set the mRGCMosaic noise flag
    theMRGCMosaic.noiseFlag = 'random';
 
    %% Set the intrinsic noise of the mRGCMosaic to 0.
    theMRGCMosaic.vMembraneGaussianNoiseSigma = 0.0;

    %% Compute noisyInstancesNum of mRGC mosaic responses with the cone mosaic noise being the only noise source
    % To do so, we compute using the input cone mosaic noisy modulation response instances
    [~, noisyResponseInstancesConeNoiseOnly] = theMRGCMosaic.compute( ...
             excitationsToModulations(theConeMosaicNoisyResponseInstances), ...
             theConeMosaicResponseTemporalSupportSeconds, ...
             'nTrials', noisyInstancesNum);

    %% Set the intrinsic noise of the mRGCMosaic to the examined level
    theMRGCMosaic.vMembraneGaussianNoiseSigma = intrinsicMRGCnoiseSigma;

    %% Compute noisyInstancesNum of mRGC mosaic responses with the mRGCMosaic intrinsic noise being the only noise source
    % To do so, we compute using the input cone mosaic noise-free modulation response 
    [~, noisyResponseInstancesIntrinsicMRGCnoiseOnly] = theMRGCMosaic.compute( ...
             excitationsToModulations(theConeMosaicNoiseFreeResponse), ...
             theConeMosaicResponseTemporalSupportSeconds, ...
             'nTrials', noisyInstancesNum);

    %% Find indices of RGCs along the y = y-stimulus position
    targetYdegs = stimulusPositionDegs(2);
    [~, mRGCIndices] = coneAndMRGCindicesAlongDesiredYposition(theMRGCMosaic, targetYdegs);

    %% Visualize response components
    generatePDFs = ~true;
    visualizeMRGMosaicResponseComponents(theMRGCMosaic, mRGCIndices, ...
        noiseFreeResponse, ...
        noisyResponseInstancesConeNoiseOnly, ...
        noisyResponseInstancesIntrinsicMRGCnoiseOnly, generatePDFs);
end