classdef fixationalEM < handle
% Fixational eye movements with realistic drift & saccade models.
%
% Syntax:
%   fixationalEM();
%
% Description:
%    Class for generating fixational eye movements using realistic drift
%    and saccade models.
%
%    For this class, The drift EM component is generated using the model
%    described in: Mergenthaler and Engbert (2007) Modeling the control of
%    fixational eye movements with neurophysiological delays", Physical
%    Review Letters, 98 (13). The microsaccde component is generated 
%    using statistical properties for the inter-saccade interval, and 
%    either a recent path avoidance strategy or a stats-based approach for
%    the micro-saccade jump.
%
% Inputs:
%   None
%
% Outputs:
%   em   -   The created fixational EM object
%
% Optional key/value pairs:
%    None
%
% See also
%   emOBJ.compute(durationSeconds, sampleTimeSeconds);
%
% References:
%    Physical Review Letters, Mergenthaler and Engbert (2007). P98.
%

% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  BW, NPC  Comments
%    06/25/18  dhb  Comment formatting.

% Public properties (user-settable)
properties
    % controlNoiseMean - Numeric.
    controlNoiseMean;

    % controlNoiseStd - Numeric.
    controlNoiseStd;

    % controlGamma - Numeric.
    controlGamma;

    % feedbackGain - Numeric.
    feedbackGain;

    % feedbackSteepness - Numeric.
    feedbackSteepness;

    % feedbackXposDelaySeconds - Numeric.
    feedbackXposDelaySeconds;

    % feedbackYposDelaySeconds - Numeric.
    feedbackYposDelaySeconds;

    % positionNoiseMean - Numeric.
    positionNoiseMean;

    % positionNoiseStd - Numeric.
    positionNoiseStd;

    % stabilizationSeconds - Numeric. How long to warm up the drift model.
    stabilizationSeconds;

    %% Microsaccade-related properties

    % microSaccadeType - String. Type of microsaccade
    % choose between {'none', 'heatmap/fixation based', 'stats based'}
    microSaccadeType;

    % microSaccadeMeanIntervalSeconds - Numeric. Mean Gamma dis of msac int
    % To clarify, the mean gamma distribution of the microsaccade intervals
    microSaccadeMeanIntervalSeconds;

    % microSaccadeMeanAmplitudeArcMin - Numeric. Mean Gamma dis of msac amp
    % To clarify, mean gamma distribution of the microsaccade amplitudes
    microSaccadeMeanAmplitudeArcMin;

    % microSaccadeMeanSpeedDegsPerSecond - Numeric. mean of norm dis speeds
    % To clarify, mean of the Normal distribution of microsaccade speeds
    microSaccadeMeanSpeedDegsPerSecond;

    % microSaccadeStDevSpeedDegsPerSecond - Numeric. STD Norm dis spds.
    % To clarify, the standard deviation of a normal distribution of
    % microsaccade speeds, in degrees per second.
    microSaccadeStDevSpeedDegsPerSecond;

    % microSaccadeTargetJitterArcMin - Numeric. MicroSaccade target jitter
    microSaccadeTargetJitterArcMin;

    % microSaccadeDirectionJitterDegs - Numeric. Directional jitter
    % The jitter in microsaccade direction, which is only used in the
    % stats-only mode.
    microSaccadeDirectionJitterDegs;

    % microSaccadeMinDurationMilliSecondsminumum duration of a microsaccade
    microSaccadeMinDurationMilliSeconds;

    %% Fixation map properties

    % fixationMapSpaceConstantArcMin - Numeric. Chance of position from 0,0
    % To expand, this is the sigma of Gaussian specifying likelihood of a
    % saccade position (from 0,0).
    fixationMapSpaceConstantArcMin;

    %% Heat map properties

    % heatMapweight - Numeric. Weight of the heat map vs the fixation map.
    heatMapWeight;

    % heatMapWidthArcMin - Numeric. The width of the heat map.
    heatMapWidthArcMin;

    % heatMapSpatialSampleArcMin - Numeric. Heat map's spatial resolution.
    heatMapSpatialSampleArcMin;

    % heatMapTemporalSampleSeconds - Numeric. Heat map temporal resolution
    heatMapTemporalSampleSeconds;

    % heatMapKernelSpaceConstantArcMin - Numeric. Spatial kernel constant.
    % To clarify, the Space constant for the heat map's spatial kernel.
    heatMapKernelSpaceConstantArcMin;

    % heatMapKernelTimeConstantSeconds - Numeric. Temporal kernel constant.
    % To clarify, the Time constant for the heat map's temporal kernel.
    heatMapKernelTimeConstantSeconds;

    % randomSeed - Numeric. The random seed to use.
    randomSeed;

    %% user-feedback

    % displayComputeProgress - Boolean. Display progress bar if set to true
    % Default is false(0). This is a numeric boolean based on observation.
    displayComputeProgress;

    % beVerbose - Boolean. Curtail output if set to false(0).
    % Default is true(1). This is a numeric boolean based on observation.
    beVerbose;
end

% Read-only properties
properties (SetAccess = private)
    % timeAxis - Vector. Time axis of the em path. Default [].
    timeAxis = [];

    % velocityArcMin - Matrix. Matrix of velocity in units of arc min/sec.
    % Matrix dimensions: nTrials x emPath length. Default [].
    velocityArcMin = [];

    % emPosArcMin - Matrix. Matrix of eye pos in arc min with no truncation
    % Matrix dimensions: nTrials x emPath length x 2. Default [].
    emPosArcMin = [];

    % emPosMicrons - Matrix. Matrix of eye positions in microns. Truncated.
    % Matrix dimensions: nTrials x emPath length x 2. Default []. 
    emPosMicrons  = [];
    
    % emPos - Matrix. Matrix of eye pos. in cone mosaic pattern size units.
    % Matrix dimensions: nTrials x emPath length x 2.
    % Note:  not set if we are computing for the new @cMosaic
    emPos = [];
end

properties (GetAccess = protected)
    % The properties contain data that are usefull for visualizing (e.g.,
    % demo app) Their values can be obtained via the subclass method
    % @visualizableFixationalEM.getValueOfProperty(propertyName)

    % emPosTimeSeriesArcMin - Last trial of emPath
    emPosTimeSeriesArcMin = [];

    % velocityArcMinPerSecTimeSeries - Last trial of Velocity
    velocityArcMinPerSecTimeSeries = [];

    % microSaccadeIntervalMilliSecondsDist - Distribution of microsac ints.
    % The distribution of the microsaccade intervals is useful for many
    % things, including visualization.
    microSaccadeIntervalsMilliSecondsDist = [];

    % microSaccadeIntervalMilliSecondsList - List with possible msac ints.
    % The list of the possible microsaccade intervals.
    microSaccadeIntervalsMilliSecondsList = [];

    % microSaccadeIntervalMilliseconds - List with generated msac ints.
    % The list of the generated microsaccade intervals.
    microSaccadeIntervalsMilliseconds = [];

    % microSaccadeAplitudesArcMinDist - Distribution of msac amplitudes.
    % The distribution of microsaccade amplitudes, useful for visualization
    microSaccadeAmplitudesArcMinDist = [];

    % microSaccadeAmplitudesArcMinList - List with possible msac amplitudes
    % The list of the possible microsaccade amplitudes.
    microSaccadeAmplitudesArcMinList = [];

    % microSaccadeAmplitudesArcMin - List with generated msac amplitudes.
    % The list of generated microsaccade amplitudes.
    microSaccadeAmplitudesArcMin = [];

    % microSaccadeSpeedDegsPerSecondDist - Dist of msac speeds in deg/sec.
    % The distribution of microsaccade speeds is useful for visualization.
    microSaccadeSpeedDegsPerSecondDist = [];

    % microSaccadeSpeedDegsPerSecondList - List of possible msac speeds.
    % The list of possible microsaccade speeds in degrees per second.
    microSaccadeSpeedDegsPerSecondList = [];

    % microSaccadeSpeedDegsPerSecond - List of generated msac speeds.
    % The list of generated microsaccade speeds in degrees per second.
    microSaccadeSpeedsDegsPerSecond = [];

    % microSaccadeDurationMillisecondsDist - Distribution of msac durations
    % The distribution of microsaccade durations (in milliseconds) is
    % derived from the speed and amplitude distributions.
    microSaccadeDurationsMillisecondsDist = [];

    % microSaccadeDurationMillisecondsList - Generated msac durations list.
    % The list of the generated microsaccade durations in milliseconds.
    microSaccadeDurationsMillisecondsList = [];

    % microSaccadeOnsetStepIndices - Generated msac onset time steps list.
    % The list of generated microsaccade onset time steps.
    microSaccadeOnsetStepIndices = [];

    % microSaccadeTargetLikelihoodMapsList - Gen. msac target chance list.
    % The list with the generated microsaccade target likelihood maps.
    microSaccadeTargetLikelihoodMapsList = [];

    % microSaccadePositions - Microsaccade positions list.
    microSaccadePositions = [];

    % heatMapTimeSeries - The heat map generated by the emPath
    heatMapTimeSeries = [];

    % heatMapXYsupport - The spatial support for the heat map.
    heatMapXYsupport = [];
end

% Properties invisible to the user
properties (Access = private)
    % stabilizationStepsNum - Numeric. The number of stabilization steps.
    stabilizationStepsNum = [];

    % tStepsNum - Numeric. Total # of time steps for entire drift duration.
    tStepsNum = [];

    % feedbackXposDelayTimeSteps - Numeric. Feedback delay for x-position.
    % Specifically, the feedback delay for the x-pos in the time steps.
    feedbackXposDelayTimeSteps;

    % feedbackYposDelayTimeSteps - Numeric. Feedback delay for y-position.
    % Specifically, the feedback delay for the y-pos in the time steps.
    feedbackYposDelayTimeSteps;

    % heatMapTimeSeriesIntermediate - Numeric. Temporary storage for the 
    % spatially-smoothed accumulated heat map of eye positions
    heatMapTimeSeriesIntermediate = [];

    % emAccumMap - Numeric. Accumulated heat map of eye positions
    emAccumMap = [];

    % heatMapUpdateIntervalStepsNum - Update interval for heat map in units
    % of time steps
    heatMapUpdateIntervalStepsNum = [];

    % lastMicroSaccadeTimeStep - Time step at which the last microsaccade 
    % occurred
    lastMicroSaccadeTimeStep = [];

    % emPosTimeSeries - The time series of emPositions
    % This is preallocated in initDriftState
    emPosTimeSeries= [];

    % velocityTimeSeries - The time series of the velocity
    velocityTimeSeries = [];

    %% Microsaccade - related

    % microSaccadeIntervalGammaShapeParameter - msac int shape param.
    % To clarify, this is the shape parameter of the Gamma distribution of
    % the microsaccade intervals.
    microSaccadeIntervalGammaShapeParameter;

    % microSaccadeAmplitudeGammaShapeParameter - msac amp shape param.
    % To clarify, this is the shape parameter of the Gamma distribution of
    % the microsaccade amplitudes.
    microSaccadeAmplitudeGammaShapeParameter;

    % microSaccadeDurationGammaShapeParameter - msac dur shape param.
    % To clarify, this is the shape parameter of the Gamma distribution of
    % the microsaccade durations.
    microSaccadeDurationGammaShapeParameter;

    % microSaccadeTargetLikelihoodSpatialMap - saccade likelihood map.
    % This is the likelihood map for saccades (Gaussian centered at (0,0))
    % with sigma : obj.fixationMapSpaceConstantArcMin
    microSaccadeTargetLikelihoodSpatialMap = [];

    % microSaccadeResidualPath - multi-step microsaccade path
    microSaccadeResidualPath = [];

    % validMicrosaccadeAmplitudeIndices - List of already-used amplitudes.
    % Reminder: sample with no replacement.
    validMicrosaccadeAmplitudeIndices = [];

    %% Drift model internal variables time series

    % positionalNoiseTimeSeries - Positional noise time series, last trial.
    % This is precomputed in initDriftState.
    positionalNoiseTimeSeries = [];

    % controlNoiseTimeSeries - Time series of control signal noise
    % This is precomputed in initDriftState.
    controlNoiseTimeSeries = [];

    % controlSignalTimeSeries - Time series of the control signal.
    % This is preallocated in initDriftState.
    controlSignalTimeSeries = [];

    % feedbackSignalTimeSeries - Time series of the feedback signal.
    % This is preallocated in initDriftState.
    feedbackSignalTimeSeries = [];

    % currentPositionHeatLevelTimeSeries - Current position heat level.
    % this could be used to trigger a saccade, but currently not used.
    currentPositionHeatLevelTimeSeries = [];

    % heatMapSpatialKernel - The heat map for the spatial kernel.
    heatMapSpatialKernel = [];

    % heatMapSpatialKernelSupport - The support for spatial kernel heat map
    heatMapSpatialKernelSupport = [];

    % heatMapTemporalKernel - The temporal kernel heat map
    heatMapTemporalKernel = [];

    % heatMapTemporalKernelSupport - Support for temporal kernel heat map.
    heatMapTemporalKernelSupport = [];
end

properties (Constant)
    % timeStepDurationSeconds - The time step duration is .001 seconds.
    timeStepDurationSeconds = 1 / 1000;

    % velocityMeasurementIntervalSeconds - Velocity interval is .041 sec.
    % This is the window for measuring velocity to 41 milliseconds, as in
    % Cherici & Rucci 2012 - Precision of sustained fixation in trained and
    % untrained observers).
    velocityMeasurementIntervalSeconds = 41/1000;

    % scalcarToArcMin - The scaling factor to convert is 0.5.
    % This 0.5 scaling factor gives drift speed in the middle of the range
    % reported byCherici et al 2012.
    scalarToArcMin = 0.5;
end

methods
    % Constructor
    function obj = fixationalEM(varargin)
        % Initialize the fixational eye movements object
        % See header commments at top of file.

        setDefaultParams(obj);
    end

    % Method to use aforementioned paper as a reference to set drift params
   setDriftParamsFromMergenthalerAndEngbert2007Paper(obj);

    % Method to set the microsaccade statistics
    setMicroSaccadeStats(obj);

    % Method to set the default parameters
    setDefaultParams(obj);

    % Method to compute the cone mosaic
    computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, varargin);

    % Method to run the computation
    compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials, ...
        computeVelocity, varargin);

    % Method to compute an emPath with duration emDurationSeconds
    % sampled with a temporal resolution of sampleDurationSeconds
    computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds);
end % Public methods

methods (Static)

    % Method to compute the fixation map
    [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = computeFixationMap(...
        timeAxis, emPaths, emPosRange, emPosDelta, varargin);

    % Method to perform the displacement analysis
    [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        performDisplacementAnalysis(emPos, timeAxisSeconds, varargin);

    % Method to generate a video of fixationalEMs with a cone mosaic
    generateEMandMosaicComboVideo(...
        fixationalEMObj, coneMosaicHexObj, varargin);
end % Static methods

methods (Access = private)

    % Method to initiate the outputs
    initOutputs(obj);

    % Method to initiate the state of the object
    initState(obj, emDurationSeconds);

    % Method to check for the microsaccade epoch
    checkForMicroSaccadeEpoch(obj,tStep);

    % Method to update the residual path
    updateMicroSaccadeResidualPath(obj);

    % Method to calculate the position of maximum likelihood
    targetPosArcMin = posOfMaxLikelihood(obj,targetLikelihoodSpatialMap);

    % Method too update the heat map
    currentPositionHeatLevel = updateHeatMap(obj, tStep);

    % Method to compute velocity
    velocityTimeSeries = computeVelocity(obj, theEMpath);

    % Method to trim, recenter, and resample the time series
    trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds);

    % Method to interpolate
    outMatrix = smartInterpolation(...
        obj,inputTimeAxis, inputMatrix, outputTimeAxis);

    % Method to compute the spatial support for the heatMap
    heatMapXYsupport = spatialSupportForHeatMap(obj);

    % Method to compute the spatial kernel for the heat map
    [kernel, kernelSupport] = spatialKernelForHeatMap(obj);

    % Method to compute the temporal kernel for the heat map
    [kernel, kernelSupport] = temporalKernelForHeatMap(obj, gridNodesNum);

    % Method to compute obj.microSaccadeTargetLikelihoodSpatialMap
    microSaccadeTargetLikelihoodSpatialMap = ...
        computeMicrosaccadeTargetLikelihoodSpatialMap(...
        obj, heatMapXYsupportArcMin);
end % private methods

end % class