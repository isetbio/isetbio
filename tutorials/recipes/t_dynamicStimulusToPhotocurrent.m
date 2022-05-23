function t_dynamicStimulusToPhotocurrent
% Illustrates how to compute photocurrent responses to a dynamic stimulus.
%
% Syntax:
%   t_dynamicStimulusToPhotocurrent
%
% Description:
%    Demonstrates how to compute photocurrent responses to a dynamic
%    stimulus whose contrast is modulated against an adapting background,
%    with realistic response noise and eye movements. If you were trying
%    to model a computational observer.
%
%    Exposes various key isetbio components and how to use them.
%
%    Key isetbio components exposed:
%      - arbitrary stimulus generated on a computer
%      - optical image sequence
%      - coneMosaic and outerSegment
%      - eyeMovements
%
%    This runs for a while and then produces three figures.
%      - Figure 1. This shows the retinal image sequence of the background
%        (uniform field) stimulus in each computed frame. The graph at the
%        upper right shows how much of a "modulation" has been combined
%        with the background, but in creating this sequence the
%        "modulation" was set equal to the background, so we just get a
%        uniform field at each timepoint. (See comments below where this
%        retinal image is created for more explanation.)  The x-axis of the
%        graph is time in seconds, and the time of each frame is given in
%        the image panels. The frames are at 60 Hz, to model a typical
%        monitor frame rate.
%
%      - Figure 2. This is the same as Figure 1, but for the modulation
%        around the background.
%
%      - Figure 3. Summarizes photocurrent responses to the modulation.
%        Top row panels display photocurrent traces separately for all L-,
%        M-, and S-cones to the adapting and the test stimulus for a single
%        trial. Note that there is no visible response modulation for the
%        test stimulus. Bottom row panels display the averaged (across
%        cones for a single trial) photocurrent responses. There are two
%        averaged responses, once across cones which respond with an
%        increased modulation to the Gabor stimulus and one across cones
%        which respond with a decreased modulation to the Gabor stimulus.
%        The noise-free averages across the same set of cones are also
%        depicted in dashed lines. These simulations demonstrate that
%        spatial response pooling across cones can be used to extract a
%        reliable estimate of the underlying response modulation even
%        though individual cone responses appear to be burried in noise.
%        The response onset latency, overshoot (150 ms) and mild adaptation
%        (200-400 ms) are due to the temporal dynamics of the photocurrent
%        impulse response. These dynamics depend on stimulus mean luminance
%        as well as cone eccentricity. In the panel titles, "adapt" refers
%        to the background alone, while "test" refers to the moduation on
%        the background. In the legend for the bottom panels, the numbers
%        in parentheses are the number of cones averaged to produce the
%        corresponding trace.
%
%    You may need to increase the size of your Matlab Java heap to run
%    this. You do that in Preferences:Matlab:Java Heap Memory, at least in
%    2017a under OS/X. Set to something bigger than the default.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/17  NPC  ISETBIO Team, 2017
%    09/01/17  npc  Wrote it.
%    09/09/17  dhb  Cosmetic pass, more comments.
%    09/12/17  npc  Superimpose eye movements on optical image sequence
%                   frames; LUT bugfix. Allow some time for respose to
%                   start returning to baseline, more comments.
%    10/02/18  jnm  Formatting

%% Initialize
ieInit;

%% Basic parameters
% eccDegs       - Mosaic eccentricity in visual degrees (> 10 results in
%                 peripheral os dynamics). 5 degrees here.
% fov           - Field of view in degrees. 0.7 degrees here.
% meanLuminance - Stimulus mean luminance (cd/m2). If you change the mean
%                 luminance, pay attention to the primary values so that
%                 they are not too close to 0 or 1. 100 Here.
% nTrials       - Response instances to compute. A single response below.
% testContrasts - Stimulus contrasts to examine. 0.2 here.
eccDegs = 5;
fov = 0.7;
meanLuminance = 100;
nTrials = 1;
testContrasts = [0.2];

%% Set the random seed
rng(1);

%% Default params for different components
% See support routines at bottom of this tutorial.
[spatialParams, temporalParams, colorModulationParams, ...
    backgroundParams] = getStimParams(fov, meanLuminance);
[oiParams, mosaicParams] = getPhysioParams(fov, eccDegs);

%% Generate a generic rectangular cone mosaic
% Start with default parameters and do some customization.
%
% center           - The mosaic's center (located on horizontal meridian),
%                    specified in meters (there are 300 microns/degree).
% integrationTime  - The mosaic's integration time
% noiseFlag        - Isomerization noise
% os               - Linear outer-segment with dynamics consistent with the
%                    mosaic's eccentricity.
% os.noiseFlag     - Photocurrent noise

fprintf('Setting up cone mosaic.\n');
cm = coneMosaic();
cm.setSizeToFOV(mosaicParams.fov);
cm.center = [eccDegs * 300, 0] * 1e-6;
cm.integrationTime = mosaicParams.integrationTimeInSeconds;
cm.noiseFlag = mosaicParams.isomerizationNoise;
cm.os = osLinear('eccentricity', eccDegs);
cm.os.noiseFlag = mosaicParams.osNoise;

%% Visualization parameters
% contrastVisualized      - Only visualize responses to the first contrast
% trialVisualized         - Only visualize the 1st response instance
% maxConesVisualized      - Visualize responses from up to this many cones
%                           for each of the L-, M- and S-cone types.
% photocurrentModlulationThresholdForInclusionToMeanResponse
%                         - In the computation of inc- and dec- response
%                           means only include cones whose modulation is >
%                           threshold x maxResponse.
% photocurrentRange       - Visualization response range for photocurrents
%                           (in pAmps)
% meanLMPhotocurrentRange - Visualization response range for modulated
%                           photocurrents from L-, and M-cones.
% meanSPhotocurrentRange  - Visualization response range for modulated
%                           photocurrents from S-cones.
contrastVisualized = 1;
trialVisualized = 1;
maxConesVisualized = 5000;
photocurrentModlulationThresholdForInclusionToMeanResponse = 0.2;
photocurrentRange = [-75 -20];
meanLMPhotocurrentRange = 1.2 * [-1 1];
meanSPhotocurrentRange = 1.2 * [-1 1];

%% Create the background scene (zero contrast)
fprintf('Creating zero contrast background scene.\n');
theBaseColorModulationParams = colorModulationParams;
theBaseColorModulationParams.coneContrasts = [0 0 0]';
theBaseColorModulationParams.contrast = 0;
backgroundScene = gaborSceneCreate(spatialParams, ...
    backgroundParams, theBaseColorModulationParams);

%% Compute the background optical image (OI).
% Optical image is the isetbio name for retinal image.
fprintf('Computing retinal image for background scene.\n');
BaseOI = opticalImageConstruct(oiParams);
oiBackground = BaseOI;
oiBackground = oiCompute(oiBackground, backgroundScene);
oiModulated = oiBackground;

%% Compute the sequence of background optical images (adapting stimulus)
% This calls a general routine that produces a blend between the background
% and a modulation. But we have set the modulation to the background here,
% so that we just get the background in every frame. The more general
% purpose of the oiSequence routine is to modulate a signal around the
% background, which will be illustrated below.
fprintf('Computing sequence of retinal images for background.\n');
zeroContrastOIsequence = oiSequence(oiBackground, oiModulated, ...
    temporalParams.sampleTimes, temporalParams.TemporalWindow, ...
    'composition', 'blend');
zeroContrastOIsequence.visualize('montage', 'showIlluminanceMap', true);

%% Generate the eye movement path for the adapting stimulus sequence
fprintf('Creating eye movement path\n');
eyeMovementsNum = ...
    zeroContrastOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
    cm.integrationTime);
theZeroContrastEMpaths = colorDetectMultiTrialEMPathGenerate(cm, 1, ...
    eyeMovementsNum, 'none');

%% Computations:
% Compute the mosaic response to the adapting stimulus, the outer segment
% IR functions and the mean photocurrents. These outer segment IR functions
% and mean photocurrents are used later on in the computation of the
% responses to the different test stimuli.
fprintf(strcat("Compute mosaic cone responses (isomerizations and ", ...
    "photocurrents) to the background sequence, taking eye ", ...
    "movements into account.\n"));
[isomerizationsAdaptingStim, photocurrentsAdaptingStim, ...
    osImpulseResponseFunctions, osMeanCurrents] = ...
    cm.computeForOISequence(zeroContrastOIsequence, ...
        'emPaths', theZeroContrastEMpaths, 'interpFilters', [], ...
        'meanCur', [], 'currentFlag', true);

%% Compute responses to all test stimuli (here different contrasts)
fprintf(strcat("Computing sequence of retinal images for background ", ...
    "plus modulated Gabor pattern at each test contrast.\n"));
fprintf(strcat("\tParameters currently set to compute for %d test ", ...
    "contrasts\n"), numel(testContrasts));
for iContrast = 1:numel(testContrasts)
    % Create the modulated scene
    colorModulationParams.contrast = testContrasts(iContrast);
    modulatedScene = gaborSceneCreate(spatialParams, ...
        backgroundParams, colorModulationParams);

    % Compute the optical image of the modulated stimulus
    oiModulated = BaseOI;
    oiModulated = oiCompute(oiModulated, modulatedScene);

    % Compute the stimulus of optical sequences
    stimulusOIsequence = oiSequence(oiBackground, oiModulated, ...
        temporalParams.sampleTimes, temporalParams.TemporalWindow, ...
        'composition', 'blend');

    % Generate the eye movement paths for all response instances
    eyeMovementsNum = ...
        stimulusOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
        cm.integrationTime);
    theEMpaths = colorDetectMultiTrialEMPathGenerate(cm, nTrials, ...
        eyeMovementsNum, temporalParams.emPathType);

    % Compute noise-free responses - this is only used for improving the
    % response visualization. Otherwise it can be skipped.
    noiseFlags = {cm.noiseFlag, cm.os.noiseFlag};
    cm.noiseFlag = 'none';
    cm.os.noiseFlag = 'none';
    fprintf(strcat("\tTest contrast %d, computing noise-free mosaic ", ...
        "responses without eye movements for visualization\n"), iContrast);
    % The emPaths below is a zero movement path for the noise-free versions
    [isomerizationsNoiseFree{iContrast}, ...
        photocurrentsNoiseFree{iContrast}, ~, ~] = ...
        cm.computeForOISequence(stimulusOIsequence, ...
            'emPaths', 0 * theEMpaths(1, :, :), ...
            'interpFilters', osImpulseResponseFunctions, ...
            'meanCur', osMeanCurrents, 'currentFlag', true);

    % Compute the response instances to the test stimulus
    cm.noiseFlag = noiseFlags{1};
    cm.os.noiseFlag = noiseFlags{2};
    fprintf(strcat("\tTest contrast %d, computing noisy mosaic ", ...
        "responses with eye movements\n"), iContrast);
    [isomerizations{iContrast}, photocurrents{iContrast}, ~, ~] = ...
        cm.computeForOISequence(stimulusOIsequence, ...
            'emPaths', theEMpaths, ...
            'interpFilters', osImpulseResponseFunctions, ...
            'meanCur', osMeanCurrents, 'currentFlag', true);

    % Visualize the stimulus oiSequence with eye movements superimposed
     eyeMovementsData = struct('show', true, ...
         'timeAxisMillisecs', cm.timeAxis * 1000, ...
         'posMicrons', squeeze(theEMpaths(trialVisualized, :, :)) ...
         * cm.patternSampleSize(1) * 1e6);
 	stimulusOIsequence.visualize('montage', 'showIlluminanceMap', true, ...
        'eyeMovementsData', eyeMovementsData);
end  % iContrast

%% Visualize the computed responses
fprintf('Creating figures for visualization ...');

% Visualize photocurrents
visualizeResponses(cm, trialVisualized, maxConesVisualized, ...
    contrastVisualized, photocurrentsNoiseFree, photocurrents, ...
    photocurrentsAdaptingStim, osMeanCurrents, photocurrentRange, ...
    meanLMPhotocurrentRange, meanSPhotocurrentRange, ...
    photocurrentModlulationThresholdForInclusionToMeanResponse);

% Wait for figure drawing
drawnow;
fprintf('done.\n');

end

% ========================================================================
% SUPPORT ROUTINES
% ========================================================================

% This one sets up parameters for a stimulus
function [spatialParams, temporalParams, colorModulationParams, ...
    backgroundParams] = getStimParams(fov, meanLuminance)
% A helper function to set up parameters for a stimulus
%
% Syntax:
%   [spatialParams, temporalParams, colorModulationParams, ...
%       backgroundParams] = getStimParams(fov, meanLuminance)
%
% Description:
%    A helper function for t_dynamicStimulusToPhotocurrent designed to set
%    up the parameters for a stimulus.
%
% Inputs:
%    fov                   - Numeric. The field of view in degrees.
%    meanLuminance         - Numeric. The mean luminance.
%
% Outputs:
%    spatialParams         - Struct. A spatial parameter structure.
%                            Contains information such as: Field of View,
%                            Gaussian FWHMD, Cycles per Degree, Spatial
%                            Phase, Orientation, Row, Col, PH, Angle, and
%                            Viewing Distance.
%    temporalParams        - Struct. A temporal parameter structure.
%                            Contains information such as: Frame Rate,
%                            Window Tau, Stimulus Duration, Stimulus
%                            Interval, EM Path Type, Sample Times, Temporal
%                            Window, and Sample Times.
%    colorModulationParams - Struct. A color modulation parameter
%                            structure. Contains information such as:
%                            Contrast, Cone Contrast, Start WI, End WI, and
%                            Delta WI.
%
% Optional key/value pairs:
%    None.
%

    % Gabor parameters
    spatialParams.fieldOfViewDegs = fov;
    spatialParams.gaussianFWHMDegs = 0.30;
    spatialParams.cyclesPerDegree = 1.75;
    spatialParams.spatialPhaseDegs = 90;
    spatialParams.orientationDegs = 0;
    spatialParams.row = 128;
    spatialParams.col = 128;
    spatialParams.ang = 0;
    spatialParams.ph = 0;
    spatialParams.viewingDistance = 0.7; % in meters

    % Color modulation parameters
    colorModulationParams.contrast = 1;
    colorModulationParams.coneContrasts = [1 1 1]';
    % Normalize the contrast.
    colorModulationParams.coneContrasts = ...
        colorModulationParams.coneContrasts ...
        ./ norm(colorModulationParams.coneContrasts);
    colorModulationParams.startWl = 380;
    colorModulationParams.endWl = 780;
    colorModulationParams.deltaWl = 4;

    % Background parameters
    lumFactor = 1.5;
    backgroundParams.backgroundxyY = [0.30 0.33 meanLuminance]';
    backgroundParams.monitorFile = 'CRT-Dell';
    backgroundParams.leakageLum = 1.0;
    backgroundParams.lumFactor = lumFactor;

    % Temporal parameters
    temporalParams.frameRate = 60;
    temporalParams.windowTauInSeconds = 0.1;
    temporalParams.stimulusDurationInSeconds = 0.8;
    temporalParams.stimulusSamplingIntervalInSeconds = ...
        1 / temporalParams.frameRate;
    % Choose emPathType from {'random', 'frozen', 'none'}
    temporalParams.emPathType = 'random';
    [temporalParams.sampleTimes, temporalParams.TemporalWindow] = ...
        flatTopGaussianWindowCreate(temporalParams);
    temporalParams.nSampleTimes = length(temporalParams.sampleTimes);
end

% This one set up parameters for the optical image and mosaic.
function [oiParams, mosaicParams] = getPhysioParams(fov, eccDegs)
% A helper function to set up parameters for an optical image and mosaic
%
% Syntax:
%   [oiParams, mosaicParams] = getPhysioParams(fov, eccDegs)
%
% Description:
%    A helper function for t_dynamicStimulusToPhotocurrent designed to set
%    up the parameters for an optical image and mosaic.
%
% Inputs:
%    fov          - Numeric. The field of view in degrees.
%    eccDegs      - Numeric. The eccentricity in degrees.
%
% Outputs:
%    oiParams     - Struct. An optical image parameter structure. Contains
%                   information such as: Off Axis, Pupil Diameter, Lens,
%                   Blur, and Field of View.
%    mosaicParams - Struct. A mosaic parameter structure. Contains
%                   information such as: OS Noise, Isomerization Noise, OS
%                   Time Step, Integration Time, Eccentricity, Aperture
%                   Blur, Inner Segment Size, LMS Ratio, Field of View, and
%                   Cone Packing.
%
% Optional key/value pairs:
%    None.
%

    % Optics parameters
    % offAxis         - Boolean. Whether to do off axis vignetting.
    % blur            - Boolean. Whether to do optical blurring.
    % lens            - Boolean. Whether to add human lens transmittance.
    % pupilDiamMm     - Numeric. The pupil diameter in degrees.
    % fieldOfViewDegs - Numeric. The field of view in degrees.
    oiParams.offAxis = false;
    oiParams.blur = false;
    oiParams.lens = false;
    oiParams.pupilDiamMm = 3;
    oiParams.fieldOfViewDegs = fov;

    % Cone mosaic parameters
    % conePacking              - String. The type of cone packing in the
    %                            mosaic. Here, rectangular.
    % fov                      - Numeric. The field of view of the mosaic
    % LMSRatio                 - Vector. The spatial densities for L-, M-,
    %                            and S-cones.
    % innerSegmentSizeMicrons  - Numeric. The inner segment diameter in
    %                            microns (square pixels).
    % apertureBlur             - Boolean. Whether or not to blur the
    %                            stimulus by the cone aperture.
    % eccentricityDegs         - Numeric. Mosaic eccentricity in degrees.
    % integrationTimeInSeconds - Numeric. Isomerizations integration time.
    % osTimeStepInSeconds      - Numeric. The time step for photocurrent
    %                            impulse response computation.
    % isomerizationNoise       - String. The isomerization Poisson noise.
    %                            choose from {'random', 'frozen', 'none'}.
    mosaicParams.osNoise = 'frozen';
    mosaicParams.conePacking = 'rect';
    mosaicParams.fov = fov;
    mosaicParams.LMSRatio = [0.62 0.31 0.07];
    mosaicParams.innerSegmentSizeMicrons = 1.4;
    mosaicParams.apertureBlur = false;
    mosaicParams.eccentricityDegs = eccDegs;
    mosaicParams.integrationTimeInSeconds = 5 / 1000;
    mosaicParams.osTimeStepInSeconds = 0.1 / 1000;
    mosaicParams.isomerizationNoise = 'frozen';
    mosaicParams.osNoise = 'frozen';
end

% Long function that produces useful visualizations for this tutorial.
function visualizeResponses(cm, trialVisualized, maxConesVisualized, ...
    contrastVisualized, photocurrentsNoiseFree, photocurrents, ...
    photocurrentsAdaptingStim, osMeanCurrents, photocurrentRange, ...
    meanLMPhotocurrentRange, meanSPhotocurrentRange, ...
    photocurrentModlulationThresholdForInclusionToMeanResponse)
% A helper function to produce helpful visualizations for the tutorial
%
% Syntax:
%   visualizeResponses(cm, trialVisualized, maxConesVisualized, ...
%       contrastVisualized, photocurrentsNoiseFree, photocurrents, ...
%       photocurrentsAdaptingStim, osMeanCurrents, photocurrentRange, ...
%       meanLMPhotocurrentRange, meanSPhotocurrentRange, ...
%       photocurrentModlulationThresholdForInclusionToMeanResponse)
%
% Description:
%    A helper function for t_dynamicStimulusToPhotocurrent designed to
%    produce useful visualizations for the tutorial.
%
% Inputs:
%    cm                        - Object. A cone mosaic object.
%    trialVisualized           - Numeric. How many trials to visualize.
%    maxConesVisualized        - Numeric. The maximum number of cones.
%    contrastVisualized        - Numeric. The contrast to visualize.
%    photocurrentsNoiseFree    - Matrix. A 4D matrix containing the
%                                noise-free photocurrents.
%    photocurrents             - Matrix. A 4D matrix containing the
%                                photocurrents.
%    photocurrentsAdaptingStim - Matrix. A 4D matrix containing the
%                                photocurrents of adapting stimuli.
%    osMeanCurrents            - Vector. The mean currents.
%    photocurrentRange         - Vector. The photocurrent range.
%    meanLMPhotocurrentRange   - Vector. The mean LM photocurrent range.
%    meanSPhotocurrentRange    - Vector. The mean S photocurrent range.
%    photocurrentModlulationThresholdForInclusionToMeanResponse
%                               - Numeric. The modulation threshold for the
%                                 photocurrent's mean response.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    currNoiseFree = RGB2XWFormat(squeeze(...
        photocurrentsNoiseFree{contrastVisualized}(...
        trialVisualized, :, :, :)));
    curr = RGB2XWFormat(squeeze(photocurrents{contrastVisualized}(...
        trialVisualized, :, :, :)));
    currAdaptingStim = RGB2XWFormat(squeeze(...
        photocurrentsAdaptingStim(1, :, :, :)));

    [peakIncLconeIndices, peakDecLconeIndices] = ...
        peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, ...
        photocurrentModlulationThresholdForInclusionToMeanResponse, 'L');
    [peakIncMconeIndices, peakDecMconeIndices] = ...
        peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, ...
        photocurrentModlulationThresholdForInclusionToMeanResponse, 'M');
    [peakIncSconeIndices, peakDecSconeIndices] = ...
        peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, ...
        photocurrentModlulationThresholdForInclusionToMeanResponse, 'S');

    incLconesNumVisualized = ...
        min([numel(peakIncLconeIndices) maxConesVisualized]);
    incMconesNumVisualized = ...
        min([numel(peakIncMconeIndices) maxConesVisualized]);
    incSconesNumVisualized = ...
        min([numel(peakIncSconeIndices) maxConesVisualized]);

    decLconesNumVisualized = ...
        min([numel(peakDecLconeIndices) maxConesVisualized]);
    decMconesNumVisualized = ...
        min([numel(peakDecMconeIndices) maxConesVisualized]);
    decSconesNumVisualized = ...
        min([numel(peakDecSconeIndices) maxConesVisualized]);

    peakIncLconeResponsesNoiseFree = ...
        currNoiseFree(peakIncLconeIndices(1:incLconesNumVisualized), :);
    peakIncMconeResponsesNoiseFree = ...
        currNoiseFree(peakIncMconeIndices(1:incMconesNumVisualized), :);
    peakIncSconeResponsesNoiseFree = ...
        currNoiseFree(peakIncSconeIndices(1:incSconesNumVisualized), :);

    peakDecLconeResponsesNoiseFree = ...
        currNoiseFree(peakDecLconeIndices(1:decLconesNumVisualized), :);
    peakDecMconeResponsesNoiseFree = ...
        currNoiseFree(peakDecMconeIndices(1:decMconesNumVisualized), :);
    peakDecSconeResponsesNoiseFree = ...
        currNoiseFree(peakDecSconeIndices(1:decSconesNumVisualized), :);

    peakIncLconeResponses = ...
        curr(peakIncLconeIndices(1:incLconesNumVisualized), :);
    peakIncMconeResponses = ...
        curr(peakIncMconeIndices(1:incMconesNumVisualized), :);
    peakIncSconeResponses = ...
        curr(peakIncSconeIndices(1:incSconesNumVisualized), :);

    peakDecLconeResponses = ...
        curr(peakDecLconeIndices(1:decLconesNumVisualized), :);
    peakDecMconeResponses = ...
        curr(peakDecMconeIndices(1:decMconesNumVisualized), :);
    peakDecSconeResponses = ...
        curr(peakDecSconeIndices(1:decSconesNumVisualized), :);

    peakIncLconeResponsesAdaptingStim = ...
        currAdaptingStim(peakIncLconeIndices(1:incLconesNumVisualized), :);
    peakIncMconeResponsesAdaptingStim = ...
        currAdaptingStim(peakIncMconeIndices(1:incMconesNumVisualized), :);
    peakIncSconeResponsesAdaptingStim = ...
        currAdaptingStim(peakIncSconeIndices(1:incSconesNumVisualized), :);

    peakDecLconeResponsesAdaptingStim = ...
        currAdaptingStim(peakDecLconeIndices(1:decLconesNumVisualized), :);
    peakDecMconeResponsesAdaptingStim = ...
        currAdaptingStim(peakDecMconeIndices(1:decMconesNumVisualized), :);
    peakDecSconeResponsesAdaptingStim = ...
        currAdaptingStim(peakDecSconeIndices(1:decSconesNumVisualized), :);

    timeAxis = cm.timeAxis;

    LincColor = [1 0.3 0.3];
    LdecColor = [0.6 0 0];

    MincColor = [0.3 0.8 0.3];
    MdecColor = [0.0 0.4 0.0];

    SincColor = [0.3 0.3 1];
    SdecColor = [0.0 0 0.6];

    subplotPosVectors = NicePlot.getSubPlotPosVectors('rowsNum', 2, ...
        'colsNum', 6, ...
        'heightMargin', 0.08, ...
        'widthMargin', 0.03, ...
        'leftMargin', 0.04, ...
        'rightMargin', 0.0, ...
        'bottomMargin', 0.06, ...
        'topMargin', 0.03);

    hFig = figure();
    clf;
    set(hFig, 'Position', [10 10 2000 900], 'Color', [1 1 1]);
    subplot('Position', subplotPosVectors(1, 1).v);
    plot(timeAxis, peakIncLconeResponsesAdaptingStim', 'r-');
    hold on;
    plot(timeAxis, peakDecLconeResponsesAdaptingStim', 'r-');
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    ylabel('photocurrent (pAmps)');
    title('L cone responses (adapt)');

    subplot('Position', subplotPosVectors(1, 2).v);
    plot(timeAxis, peakIncLconeResponses', 'r-');
    hold on;
    plot(timeAxis, peakDecLconeResponses', 'r-');
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    title('L cone responses (test)');

    subplot('Position', subplotPosVectors(1, 3).v);
    plot(timeAxis, peakIncMconeResponsesAdaptingStim', '-', ...
        'Color', [0 0.8 0.0]);
    hold on;
    plot(timeAxis, peakDecMconeResponsesAdaptingStim', '-', ...
        'Color', [0 0.8 0.0]);
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    title('M cone responses (adapt)');

    subplot('Position', subplotPosVectors(1, 4).v);
    plot(timeAxis, peakIncMconeResponses', '-', 'Color', [0 0.7 0.0]);
    hold on;
    plot(timeAxis, peakDecMconeResponses', '-', 'Color', [0 0.7 0.0]);
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    title('M cone responses (test)');

    subplot('Position', subplotPosVectors(1, 5).v);
    plot(timeAxis, peakIncSconeResponsesAdaptingStim', 'b-');
    hold on;
    plot(timeAxis, peakDecSconeResponsesAdaptingStim', 'b-');
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    title('S cone responses (adapt)');

    subplot('Position', subplotPosVectors(1, 6).v);
    plot(timeAxis, peakIncSconeResponses', 'b-');
    hold on;
    plot(timeAxis, peakDecSconeResponses', 'b-');
    set(gca, 'YLim', photocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:5:0, 'FontSize', 14);
    grid on;
    box off;
    title('S cone responses (test)');

    subplot('Position', subplotPosVectors(2, 1).v);
    plot(timeAxis, mean(peakIncLconeResponsesAdaptingStim, 1), '-', ...
        'Color', LincColor, 'LineWidth', 1.5);
    hold on;
    plot(timeAxis, mean(peakDecLconeResponsesAdaptingStim, 1), '-', ...
        'Color', LdecColor, 'LineWidth', 1.5);
    set(gca, 'YLim', osMeanCurrents(1) + meanLMPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    grid on;
    box off;
    ylabel('photocurrent (pAmps)');
    xlabel('time (seconds)');
    title('mean L cone response (adapt)');

    subplot('Position', subplotPosVectors(2, 2).v);
    plot(timeAxis, mean(peakIncLconeResponsesNoiseFree, 1), '-', ...
        'Color', LincColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    hold on
    plot(timeAxis, mean(peakDecLconeResponsesNoiseFree, 1), '-', ...
        'Color', LdecColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    plot(timeAxis, mean(peakIncLconeResponses, 1), '-', ...
        'Color', LincColor, 'LineWidth', 1.5);
    plot(timeAxis, mean(peakDecLconeResponses, 1), '-', ...
        'Color', LdecColor, 'LineWidth', 1.5);

    set(gca, 'YLim', osMeanCurrents(1) + meanLMPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    legend({'mean-inc, no noise', 'mean-dec, no noise', ...
        sprintf('mean-inc (%d)', incLconesNumVisualized), ...
        sprintf('mean-dec (%d)', decLconesNumVisualized)}, ...
        'Location', 'SouthEast');
    grid on;
    box off;
    title('mean L cone response (test)');

    subplot('Position', subplotPosVectors(2, 3).v);
    plot(timeAxis, mean(peakIncMconeResponsesAdaptingStim, 1), '-', ...
        'Color', MincColor, 'LineWidth', 1.5);
    hold on;
    plot(timeAxis, mean(peakDecMconeResponsesAdaptingStim, 1), '-', ...
        'Color', MdecColor, 'LineWidth', 1.5);
    set(gca, 'YLim', osMeanCurrents(2) + meanLMPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    grid on;
    box off;
    title('mean M cone response (adapt)');

    subplot('Position', subplotPosVectors(2, 4).v);
    plot(timeAxis, mean(peakIncMconeResponsesNoiseFree, 1), '-', ...
        'Color', MincColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    hold on;
    plot(timeAxis, mean(peakDecMconeResponsesNoiseFree, 1), '-', ...
        'Color', MdecColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    plot(timeAxis, mean(peakIncMconeResponses, 1), '-', ...
        'Color', MincColor, 'LineWidth', 1.5);
    plot(timeAxis, mean(peakDecMconeResponses, 1), '-', ...
        'Color', MdecColor, 'LineWidth', 1.5);
    set(gca, 'YLim', osMeanCurrents(2) + meanLMPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    legend({'mean-inc, no noise', 'mean-dec, no noise', ...
        sprintf('mean-inc (%d)', incMconesNumVisualized), ...
        sprintf('mean-dec (%d)', decMconesNumVisualized)}, ...
        'Location', 'SouthEast');

    grid on;
    box off;
    title('mean M cone response (test)');

    subplot('Position', subplotPosVectors(2, 5).v);
    plot(timeAxis, mean(peakIncSconeResponsesAdaptingStim, 1), '-', ...
        'Color', SincColor, 'LineWidth', 1.5);
    hold on;
    plot(timeAxis, mean(peakDecSconeResponsesAdaptingStim, 1), '-', ...
        'Color', SdecColor, 'LineWidth', 1.5);
    set(gca, 'YLim', osMeanCurrents(3) + meanSPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    grid on;
    box off;
    title('mean S cone response (adapt)');

    subplot('Position', subplotPosVectors(2, 6).v);
    plot(timeAxis, mean(peakIncSconeResponsesNoiseFree, 1), '-', ...
        'Color', SincColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    hold on
    plot(timeAxis, mean(peakDecSconeResponsesNoiseFree, 1), '-', ...
        'Color', SdecColor * 0.4 + 0.6 * [1 1 1], 'LineWidth', 8.0);
    plot(timeAxis, mean(peakIncSconeResponses, 1), '-', ...
        'Color', SincColor, 'LineWidth', 1.5);
    plot(timeAxis, mean(peakDecSconeResponses, 1), '-', ...
        'Color', SdecColor, 'LineWidth', 1.5);
    set(gca, 'YLim', osMeanCurrents(3) + meanSPhotocurrentRange, ...
        'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, ...
        'YTick', -100:1:0, 'FontSize', 14);
    legend({'mean-inc, no noise', 'mean-dec, no noise', ...
        sprintf('mean-inc (%d)', incSconesNumVisualized), ...
        sprintf('mean-dec (%d)', decSconesNumVisualized)}, ...
        'Location', 'SouthEast');
    grid on;
    box off;
    title('mean S cone response (test)');
end

% Function to sort cones according to the magnitude of their photocurrent
% modulations. Currently, this only workds for rect coneMosaics.
function [peakIncrementConeIndices, peakDecrementConeIndices] = ...
    peakResponseConeIndices(cm, currents, osMeanCurrents, ...
    photocurrentModlulationThresholdForInclusionToMeanResponse, coneType)
% Sort rect ConeMosaic cones by  magnitude of photocurrent modulations
%
% Syntax:
%   [peakIncrementConeIndices, peakDecrementConeIndices] = ...
%       peakResponseConeIndices(cm, currents, osMeanCurrents, ...
%       photocurrentModlulationThresholdForInclusionToMeanResponse, ...
%       coneType)
%
% Description:
%    This is a helper function in t_dynamicStimulusToPhotocurrent, designed
%    to sort cones according to the magnitude of their photocurrent
%    modulations. At present it only works for rectangular coneMosaics.
%
% Inputs:
%    cm                       - Object. A cone mosaic object.
%    currents                 - Matrix. The currents.
%    osMeanCurrents           - Vector. The L, M, and S mean currents.
%    photocurrentModlulationThresholdForInclusionToMeanResponse
%                             - Numeric. The modulation threshold for the
%                               photocurrent to serve as the limit for
%                               inclusion into the mean response.
%    coneType                 - Char. The character indicating the cone
%                               type. Options are 'L', 'M', and 'S'.
%
% Outputs:
%    peakIncrementConeIndices - Vector. The peak increment cone indices.
%    peakDecrementConeIndices - Vector. The peak decrement cone indices.
%
% Optional key/value pairs:
%    None.
%
    if isa(cm, 'coneMosaicHex')
        error('This method only works for rectangular cone mosaics');
    end
    % Determine target cone index
    [~, coneIndex] = ismember(coneType, {'L', 'M', 'S'});
    if coneIndex == 0
        peakIncrementConeIndices = [];
        peakDecrementConeIndices = [];
        return;
    end

    % Compute response modulations for the target cone type
    pattern = cm.pattern(:);
    targetConeIndices = find(pattern == coneIndex + 1);
    currentModulations = currents(targetConeIndices, :) ...
        - osMeanCurrents(coneIndex);
    thresholdModulation = ...
        photocurrentModlulationThresholdForInclusionToMeanResponse ...
        * max(abs(currentModulations(:)));

    % Sort cones with positive modulation
    m = max(currentModulations, [], 2);
    [m, idx] = sort(m, 'descend');
    idx = idx(m > thresholdModulation);
    peakIncrementConeIndices = targetConeIndices(idx);

    % Sort cones with negative modulation
    m = min(currentModulations, [], 2);
    [m, idx] = sort(m, 'ascend');
    idx = idx(m < -thresholdModulation);
    peakDecrementConeIndices = targetConeIndices(idx);
end

% Function to generate a flat-top Gaussian temporal modulation
function [sampleTimes, TemporalWindow] = ...
    flatTopGaussianWindowCreate(temporalParams)
% Helper function to generate a flat-top gaussian temporal modulation
%
% Syntax:
%   [sampleTimes, TemporalWindow] = ...
%       flatTopGaussianWindowCreate(temporalParams)
%
% Description:
%    This is a helper function for t_dynamicStimulusToPhotocurrent designed
%    to generate a flat-top gaussian temporal modulation.
%
% Inputs:
%    temporalParams - Struct. A structure containing the necessary temporal
%                     parameters for creating the window. Contains:
%           stimulusSamplingIntervalInSeconds: Numeric. See name.
%           stimulusDurationInSeconds: Numeric. See name.
%
% Outputs:
%    sampleTimes    - Vector. The sample times in a row vector.
%    TemporalWindow - Vector. The column vector to match sample times.
%
% Optional key/value pairs:
%    None.
%
    sampleTimes = 0:temporalParams.stimulusSamplingIntervalInSeconds:...
        temporalParams.stimulusDurationInSeconds;
    returnToBaseTimeSeconds = 0.1;
    returnToBaseLineSamples = round(returnToBaseTimeSeconds / ...
        temporalParams.stimulusSamplingIntervalInSeconds);
    L = numel(sampleTimes);
    upRampTime = temporalParams.windowTauInSeconds;
    upRampSamples = ceil(upRampTime / ...
        temporalParams.stimulusSamplingIntervalInSeconds);
    downRampSamples = L - upRampSamples + 1 - returnToBaseLineSamples;
    TemporalWindow = ones(L, 1);
    TemporalWindow(1:upRampSamples) = 1 - exp(-0.5 * ...
        (linspace(0, 1, upRampSamples) / 0.33) .^ 2);
    TemporalWindow(downRampSamples - returnToBaseLineSamples + ...
        (1:upRampSamples)) = TemporalWindow(upRampSamples:-1:1);
    TemporalWindow(downRampSamples - returnToBaseLineSamples + ...
        upRampSamples + 1:L) = 0;
    visualizeWindow = false;
    if (visualizeWindow)
        figure(2);
        plot(sampleTimes, TemporalWindow, 'ks-');
        drawnow;
    end
end

% Function to compute the XYZtoConeExcitation matrix
function [M, coneFundamentals, coneSpectralSampling] = XYZToCones
% Use Stockman-Sharpe 2deg fundamentals & CIE corresponding xyz functions
%
% Syntax:
%   [M, coneFundamentals, coneSpectralSampling] = XYZToCones
%
% Description:
%    Here we'll use the Stockman-Sharpe 2-degree fundamentals and the
%    proposed CIE corresponding XYZ functions.
%
% Inputs:
%    None.
%
% Outputs:
%    M                    - Matrix. Matrix used to help calculate the
%                           background excitations.
%    coneFundamentals     - Matrix. A matrix of cone fundamentals.
%    coneSpectralSampling - Vector. The RGB cone spectral sampling.
%
% Optional key/value pairs:
%    None.
%
    theCones = load('T_cones_ss2');
    theXYZ = load('T_xyzCIEPhys2');

    XYZcolorMatchingFunctions = 683 * theXYZ.T_xyzCIEPhys2;
    coneFundamentals = 683 * theCones.T_cones_ss2;
    coneSpectralSampling = theCones.S_cones_ss2;
    M = ((XYZcolorMatchingFunctions') \ (coneFundamentals'))';
end

% Generate a color gabor ISETBio scene with known properties.
function theScene = gaborSceneCreate(...
    spatialParams, backgroundParams, colorModulationParams)
% Generate a color gabor ISETBio scene with known properties.
%
% Syntax:
%   theScene = gaborSceneCreate(...
%       spatialParams, backgroundParams, colorModulationParams)
%
% Description:
%    Function to generate a color Gabor isetbio scene generated on a
%    display. This generates a color Gabor that has known properties with
%    respect to the Stockman-Sharpe cone fundamentals. This would be a
%    common thing to construct in a psychophysical experiment that measured
%    color detection thresholds.
%
% Inputs:
%    spatialParams         - Struct. A structure containing the following:
%           fieldOfViewDegs: Numeric. The field of view, in degrees.
%           gaussianFWHMDegs: Numeric. The full width half max Gaussian.
%           cyclesPerDegree: Numeric. The cycles per degree.
%           spatialPhaseDegs: Numeric. The spatial phase, in degrees.
%           orientationDegs: Numeric. The orientation, in degrees.
%           row: Numeric. The row.
%           col: Numeric. The column
%           ang: Numeric. The angle in degrees.
%           ph: Numeric. The phase.
%           viewingDistance: Numeric. The viewing distance, in meters.
%    backgroundParams      - Struct. A structure containing the following:
%           backgroundxyY: Vector. A 3x1 numeric column vector.
%           monitorFile: String. A string containing the monitor filename.
%           leakageLum: Numeric. The numeric value of luminance leakage.
%           lumFactor: Numeric. The luminance factor.
%    colorModulationParams - Struct. A structure containing the following:
%           contrast: Numeric. The contrast value.
%           coneContrasts: Vector. A 3x1 numeric column vector.
%           startWl: Numeric. The starting wavelength.
%           endWl: Numeric. The ending wavelength.
%           deltaWl: Numeric.  The delta for wavelengths.
%
% Outputs:
%    theScene              - Struct. A scene structure.
%
% Optional key/value pairs:
%    None.
%
    if (~isfield(spatialParams, 'fieldOfViewDegs'))
        error('Spatial parameters must have a fieldOfViewDegs field');
    end
    fieldOfViewDegs = spatialParams.fieldOfViewDegs;

    % Make the spatial pattern, and make it as a Gabor and convert to a
    % modulation around the mean (mean == 0)
    spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(...
        spatialParams, colorModulationParams.contrast));
    spatialModulation = spatialPattern - 1;

    % Make sure that contrast & background vectors are both column vectors.
    coneContrast = colorModulationParams.coneContrasts(:);
    backgroundxyY = backgroundParams.backgroundxyY(:);
    %backgroundxyY(3) = backgroundxyY(3) / backgroundParams.lumFactor;

    % Convert pattern to a color modulation specified in cone space
    [M_XYZToCones, T_cones, S_cones] = XYZToCones();

    % Convert background to cone excitations
    backgroundConeExcitations = M_XYZToCones * xyYToXYZ(backgroundxyY);

    % Convert test cone contrasts to cone excitations
    testConeExcitationsDir = (coneContrast .* backgroundConeExcitations);

    % Make the color pattern in LMS excitations
    patternConeExcitationsBg = ones(spatialParams.row, spatialParams.col);
    patternConeExcitations = ...
        zeros(spatialParams.row, spatialParams.col, 3);
    for ii = 1:3
        patternConeExcitations(:, :, ii) = patternConeExcitationsBg ...
            * backgroundConeExcitations(ii) ...
            + spatialModulation * testConeExcitationsDir(ii);
    end

    % We need a display. We'll just use the description of a CRT that we
    % have handy. In doing so, we are assuming that the differences between
    % CRT's used in different threshold experiments do not have a
    % substantial effect on the thresholds. There will be a little effect
    % because differences in channel spectra will lead to differences in
    % the retinal image because of chromatic aberration, but given the
    % general similarity of monitor channel spectra we expect these
    % differences to be small.

    display = displayCreate(backgroundParams.monitorFile);

    % Specify the SPDs
    display = displaySet(display, 'spd', ...
        backgroundParams.lumFactor * displayGet(display, 'spd'));

    % Specify a linear 16-bit LUT
    displayBits = 16;
    display = displaySet(display, 'gTable', ...
        repmat(linspace(0, 1, 2 ^ displayBits)', [1 3]));

    % Set the viewing distance
    display = displaySet(display, ...
        'viewingdistance', spatialParams.viewingDistance);

    % Reset the display wavelengths, since everything else follows on this.
    newWls = colorModulationParams.startWl: ...
        colorModulationParams.deltaWl:colorModulationParams.endWl;
    display = displaySet(display, 'wave', newWls);

    % Get display channel spectra. The S vector displayChannelS is PTB
    % format for specifying wavelength sampling:
    %    [startWl deltaWl nWlSamples]
    displayChannelWavelengths = displayGet(display, 'wave');
    displayChannelS = WlsToS(displayChannelWavelengths);
    displayChannelSpectra = displayGet(display, 'spd');

    % Spline cones to same wavelength sampling as display
    T_conesForDisplay = ...
        SplineCmf(S_cones, T_cones, displayChannelWavelengths);

    % Find the matrix that converts between linear channel weights (called
    % "primary" in PTB lingo) and LMS excitations, and its inverse.
    % Multiplication by the deltaWl is to handle fact that in isetbio
    % radiance is specified in Watts/[m2-sr-nm].
    %
    % Also get matrices for going in and out of XYZ, and compute display
    % max luminance as a sanity check.
    M_PrimaryToConeExcitations = T_conesForDisplay ...
        * displayChannelSpectra * displayChannelS(2);
    M_ConeExcitationsToPrimary = inv(M_PrimaryToConeExcitations);

    % Convert the cone excitations image to RGB
    [patternConeExcitationsCalFormat, m, n] = ...
        ImageToCalFormat(patternConeExcitations);
    patternPrimaryCalFormat = M_ConeExcitationsToPrimary ...
        * patternConeExcitationsCalFormat;
    patternPrimary = CalFormatToImage(patternPrimaryCalFormat, m, n);

    % Check that the image is within the monitor gamut. If the pattern
    % represents an actual stimulus produced with an actual monitor, things
    % should be OK if both are represented properly in this routine.
    maxPrimary = max(patternPrimaryCalFormat(:));
    minPrimary = min(patternPrimaryCalFormat(:));
    if ((minPrimary < 0.1) || (maxPrimary > 0.9))
        fprintf(2, strcat("You may want to adjust the mean luminance ", ...
            "or the lumFactor\n"));
        fprintf(2, strcat("Primary value range for the R channel: ", ...
            "%2.4f - %2.4f\n"), min(patternPrimaryCalFormat(1, :)), ...
            max(patternPrimaryCalFormat(1, :)));
        fprintf(2, strcat("Primary value range for the G channel: ", ...
            "%2.4f - %2.4f\n"), min(patternPrimaryCalFormat(2, :)), ...
            max(patternPrimaryCalFormat(2, :)));
        fprintf(2, strcat("Primary value range for the B channel: ", ...
            "%2.4f - %2.4f\n"), min(patternPrimaryCalFormat(3, :)), ...
            max(patternPrimaryCalFormat(3, :)));
    end

    if (maxPrimary > 1 || minPrimary < 0)
        error(strcat("RGB primary image is out of gamut. You need to ", ...
            "do something about this"));
    end

    % Gamma correct the primary values, so we can pop them into an isetbio
    % scene in some straightforward manner. It's important to have a lot of
    % steps in the inverse gamma, so that one doesn't truncate very low
    % contrast scenes.
    patternRGB = round(ieLUTLinear(patternPrimary, ...
        displayGet(display, 'inverse gamma', 2 ^ displayBits)));

    % Finally, make the actual isetbio scene
    % This combines the image we build and the display properties.
    theScene = sceneFromFile(patternRGB, 'rgb', [], display);
    theScene = sceneSet(theScene, 'h fov', fieldOfViewDegs);
    fprintf('Mean scene luminance: %2.2f cd/m2\n', ...
        sceneGet(theScene, 'mean lum'));
end

% Convert stimulus params to a format expected by isetbios's imageHarmonic.
function imageHarmonicParams = imageHarmonicParamsFromGaborParams(...
    spatialParams, contrast)
% Convert stimulus params to a format expected by isetbios's imageHarmonic.
%
% Syntax:
%   imageHarmonicParams = imageHarmonicParamsFromGaborParams(...
%       spatialParams, contrast)
%
% Description:
%    This function is intended to convert stimulus params to a format
%    expected by isetbios's imageHarmonic.
%
% Inputs:
%    spatialParams       - Struct. A structure containing the following:
%           fieldOfViewDegs: Numeric. The field of view, in degrees.
%           gaussianFWHMDegs: Numeric. The full width half max Gaussian.
%           cyclesPerDegree: Numeric. The cycles per degree.
%           spatialPhaseDegs: Numeric. The spatial phase, in degrees.
%           orientationDegs: Numeric. The orientation, in degrees.
%           row: Numeric. The row.
%           col: Numeric. The column
%           ang: Numeric. The angle in degrees.
%           ph: Numeric. The phase.
%           viewingDistance: Numeric. The viewing distance, in meters.
%    contrast            - Numeric. Number representing the contrast level.
%
% Outputs:
%    imageHarmonicParams - Struct. A structure of the harmonic parameters.
%                          The possible parameters are as follows:
%           ang: Orientation angle of grating. Default is 0.
%           contrast: Image contrast. Default is 1.
%           freq: Spatial frequency. Default is 1.
%           ph: Phase. Default is pi / 2
%           row: Window rows. Default is 64.
%           col: Window columns. Default is 64.
%           GaborFlag: Gaussian window, standard deviation re: window size.
%
% Optional key/value pairs:
%    None.
%

    % Set up base parameters
    imageHarmonicParams = spatialParams;

    % Pull contrast into the structure
    imageHarmonicParams.contrast = contrast;

    % Computed parameters. These convert numbers to a form used by
    % underlying routines. This one is spatial frequency
    cyclesPerImage = spatialParams.fieldOfViewDegs ...
        * spatialParams.cyclesPerDegree;
    imageHarmonicParams.freq = cyclesPerImage;

    if (isfield(spatialParams, 'spatialPhaseDegs'))
        imageHarmonicParams.ph = spatialParams.spatialPhaseDegs / 180 * pi;
    end

    if (isfield(spatialParams, 'orientationDegs'))
        imageHarmonicParams.ang = spatialParams.orientationDegs / 180 * pi;
    end

    % Set GaborFlag to specify window.
    gaussianStdDegs = FWHMToStd(spatialParams.gaussianFWHMDegs);
    gaussianStdImageFraction = gaussianStdDegs ...
        / spatialParams.fieldOfViewDegs;
    imageHarmonicParams.GaborFlag = gaussianStdImageFraction;
end

% Function to compute an optical image object
function theOI = opticalImageConstruct(oiParams)
% A function to compute an optical image
%
% Syntax:
%   theOI = opticalImageConstruct(oiParams)
%
% Description:
%    A function to compute an optical image object(structure).
%
% Inputs:
%    oiParams - Struct. A structure containing the OI Parameters, such as:
%
% Outputs:
%    theOI    - Struct. The created optical image structure. Contains:
%           type: String. Always "opticalImage"
%           name: String. A string containing the OI's name
%           data: Struct. A structure containing data (bitdepth: Numeric)
%           optics: Struct. An optical image structure.
%           diffuser: Struct. A structure containing diffuser parameters.
%           consistency: Boolean. Is the data consistent?
%
% Optional key/value pairs:
%    None.
%

    % Basic create.
    theOI = oiCreate('wvf human', oiParams.pupilDiamMm);
    theOI = oiSet(theOI, 'h fov', oiParams.fieldOfViewDegs);

    % Set the pupil diamter
    focalLength = oiGet(theOI, 'distance');
    desiredFNumber = focalLength / (oiParams.pupilDiamMm / 1000);
    theOI  = oiSet(theOI , 'optics fnumber', desiredFNumber);
    pupilDiamMmCheck = 1000 * oiGet(theOI, 'optics aperture diameter');
    if (max(abs(pupilDiamMmCheck - oiParams.pupilDiamMm)) > 1e-8)
        error('Failed to set pupil diameter as expected');
    end

    % Take out off axis vignetting if requested
    optics = oiGet(theOI, 'optics');
    if (~oiParams.offAxis)
        optics = opticsSet(optics, 'off axis method', 'skip');
    end
    theOI = oiSet(theOI, 'optics', optics);

    % Take out optical blurring if requested
    optics = oiGet(theOI, 'optics');
    if (~oiParams.blur)
        optics = opticsSet(optics, ...
            'OTF', ones(size(opticsGet(optics, 'OTF'))));
    end
    theOI = oiSet(theOI, 'optics', optics);

    % Take out lens transmittance if desired
    if (~oiParams.lens)
        lens = oiGet(theOI, 'lens');
        lens.density = 0;
        theOI = oiSet(theOI, 'lens', lens);
    end
end

% Function to create an array of EMpaths, one for each of the nTrials
function theEMpaths = colorDetectMultiTrialEMPathGenerate(...
    theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
% Function to create an array of EMpaths, one for each of the nTrials
%
% Syntax:
%   theEMpaths = colorDetectMultiTrialEMPathGenerate(theConeMosaic, ...
%       nTrials, eyeMovementsPerTrial, emPathType, [varargin])
%
% Description:
%    Function to create an array of EMpaths, one for each of the nTrials.
%
% Inputs:
%    theConeMosaic        - Object. A cone mosaic object.
%    nTrials              - Numeric. The number of trials
%    eyeMovementsPerTrial - Numeric. The number of eye movements per trial.
%    emPathType           - String. The eye movement path type.
%    varargin             - (Optional) Varies. Additional arguments, as
%                           desired, to provide information.
%
% Outputs:
%    theEMpaths           - Matrix. Two matrices, one for each eye,
%                           denoting the eye movement paths.
%
% Optional key/value pairs:
%    None.
%

    p = inputParser;
    p.addParameter('seed', 1, @isnumeric);
    p.addParameter('centeredEMPaths', false, @islogical);
    p.parse(varargin{:});

    switch (emPathType)
        case 'none'
            theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
        case 'frozen'
            rng(p.Results.seed);
            theFixedEMpath = ...
                theConeMosaic.emGenSequence(eyeMovementsPerTrial);
            theEMpaths = repmat(theFixedEMpath, [nTrials 1 1]);
        case 'frozen0'
            rng(p.Results.seed);
            theFixedEMpath = 0 ...
                * theConeMosaic.emGenSequence(eyeMovementsPerTrial);
            theEMpaths = repmat(theFixedEMpath, [nTrials 1 1]);
        case 'random'
            theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
            for iTrial = 1:nTrials
                theEMpaths(iTrial, :, :) = ...
                    theConeMosaic.emGenSequence(eyeMovementsPerTrial);
            end
        otherwise
            error(strcat("Unknown emPathType: '%s'. Valid choices: ", ...
                "'none', 'frozen', 'random'."), emPathType);
    end

    if (p.Results.centeredEMPaths)
        % find centers
        emPathCenters = mean(theEMpaths, 2);
        % subtract them from the emPaths
        theEMpaths = bsxfun(@minus, theEMpaths, emPathCenters);
        % emPaths must be integered-valued, so round
        theEMpaths = round(theEMpaths);
    end
end