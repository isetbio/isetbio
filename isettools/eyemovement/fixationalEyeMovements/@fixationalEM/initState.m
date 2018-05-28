function initState(obj, emDurationSeconds)
% Initialize a fixationalEM object
%
% Syntax:
%   obj.initState(emDurationSeconds)
%   initState(obj, emDurationSeconds)
%
% Description:
%    Initialize the state of a fixationalEM object.
%
% Inputs:
%    obj               - Object. A fixationalEM object.
%    emDurationSeconds - Numeric. The eye movement duration, in seconds.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting

% Compute stabilization steps
obj.stabilizationStepsNum = round(obj.stabilizationSeconds / ...
    obj.timeStepDurationSeconds);

% Compute total number of time steps
obj.tStepsNum = round(emDurationSeconds / obj.timeStepDurationSeconds) ...
    + obj.stabilizationStepsNum + 1;

% Pre-allocate memory for time axis          
obj.timeAxis = zeros(1, obj.tStepsNum);

% Drift-related params
% Compute feedback delay in time steps for X/Y position
obj.feedbackXposDelayTimeSteps = round(obj.feedbackXposDelaySeconds / ...
    obj.timeStepDurationSeconds);
obj.feedbackYposDelayTimeSteps = round(obj.feedbackYposDelaySeconds / ...
    obj.timeStepDurationSeconds);

% Compute positional noise for the entire duration of the drift
obj.positionalNoiseTimeSeries = normrnd(...
    obj.positionNoiseMean * ones(2, obj.tStepsNum), ...
    obj.positionNoiseStd  * ones(2, obj.tStepsNum));

% Compute control signal noise for the entire duration of the drift
obj.controlNoiseTimeSeries = normrnd(...
    obj.controlNoiseMean * ones(2, obj.tStepsNum), ...
    obj.controlNoiseStd  * ones(2, obj.tStepsNum));

% Pre-allocate memory for the feedback and the control components of the
% drift mechanism
obj.feedbackSignalTimeSeries = zeros(2, obj.tStepsNum);
obj.controlSignalTimeSeries = zeros(2, obj.tStepsNum);

% Pre-allocate memory for the computed emPositions - assuming 1 millisecond
% time sample
obj.emPosTimeSeries = zeros(2, obj.tStepsNum);

% Pre-allocate memory for velocity
obj.velocityArcMinPerSecTimeSeries = zeros(1, obj.tStepsNum);
obj.velocityTimeSeries = zeros(1, obj.tStepsNum);

% Pre-allocate memory for the computed emPositions - arc min
obj.emPosTimeSeriesArcMin = zeros(2, obj.tStepsNum);

% Time steps where a saccade is initiated
obj.microSaccadeOnsetStepIndices = [];
obj.microSaccadeAmplitudesArcMin = [];
obj.microSaccadeSpeedsDegsPerSecond = [];
obj.microSaccadeIntervalsMilliseconds = [];
obj.microSaccadeTargetLikelihoodMapsList = [];
obj.microSaccadePositions = [];

obj.heatMapXYsupport = obj.spatialSupportForHeatMap();
obj.microSaccadeTargetLikelihoodSpatialMap = ...
    obj.computeMicrosaccadeTargetLikelihoodSpatialMap(...
    obj.heatMapXYsupport);

tBins = numel(obj.stabilizationStepsNum:obj.tStepsNum);
obj.heatMapTimeSeries = zeros(tBins, numel(obj.heatMapXYsupport), ...
    numel(obj.heatMapXYsupport));
obj.currentPositionHeatLevelTimeSeries = zeros(1, tBins);

obj.heatMapUpdateIntervalStepsNum = round(...
    obj.heatMapTemporalSampleSeconds / obj.timeStepDurationSeconds);
tBinsSubSampled = 1:obj.heatMapUpdateIntervalStepsNum:...
    numel(obj.stabilizationStepsNum);
obj.heatMapTimeSeriesIntermediate = zeros(numel(tBinsSubSampled), ...
    numel(obj.heatMapXYsupport), numel(obj.heatMapXYsupport));

% Compute spatial kernel for heat map
[obj.heatMapSpatialKernel, obj.heatMapSpatialKernelSupport] = ...
    obj.spatialKernelForHeatMap();

% Compute temporal kernel for heat map
[obj.heatMapTemporalKernel, obj.heatMapTemporalKernelSupport] = ...
    obj.temporalKernelForHeatMap(length(obj.heatMapXYsupport));

% Generate Gamma distibuted microsaccade intervals
maxMicroSaccadeNum = obj.tStepsNum - obj.stabilizationStepsNum;
shape = obj.microSaccadeIntervalGammaShapeParameter;
scale = obj.microSaccadeMeanIntervalSeconds * 1000 / shape;
obj.microSaccadeIntervalsMilliSecondsDist = round(random('Gamma', ...
    shape, scale, [1 max([10000 obj.tStepsNum])]));
obj.microSaccadeIntervalsMilliSecondsList = ...
    obj.microSaccadeIntervalsMilliSecondsDist(1:maxMicroSaccadeNum);

% Generate microsaccade speed distribution
obj.microSaccadeSpeedDegsPerSecondDist = random('Normal', ...
    obj.microSaccadeMeanSpeedDegsPerSecond, ...
    obj.microSaccadeStDevSpeedDegsPerSecond, ...
    [1 max([10000 obj.tStepsNum])]);
obj.microSaccadeSpeedDegsPerSecondList = ...
    obj.microSaccadeSpeedDegsPerSecondDist(1:maxMicroSaccadeNum);

if (strcmp(obj.microSaccadeType, 'stats based'))
    % Generate Gamma distributed microsaccade amplitudes
    shape = obj.microSaccadeAmplitudeGammaShapeParameter;
    scale = obj.microSaccadeMeanAmplitudeArcMin / shape;
    obj.microSaccadeAmplitudesArcMinDist = random('Gamma', shape, ...
        scale, [1 max([10000 obj.tStepsNum])]);
    obj.microSaccadeAmplitudesArcMinList = ...
        obj.microSaccadeAmplitudesArcMinDist(1:maxMicroSaccadeNum);
    obj.validMicrosaccadeAmplitudeIndices = ...
        1:numel(obj.microSaccadeAmplitudesArcMinList);

    % Generate microsaccade durations distribution from the speed and the
    % amplitude distributions
    obj.microSaccadeDurationsMillisecondsDist = round(...
        obj.microSaccadeAmplitudesArcMinDist ./ ...
        (obj.microSaccadeSpeedDegsPerSecondDist / 1000 * 60));
    obj.microSaccadeDurationsMillisecondsDist(...
        obj.microSaccadeDurationsMillisecondsDist < ...
        obj.microSaccadeMinDurationMilliSeconds) = ...
        obj.microSaccadeMinDurationMilliSeconds;
    obj.microSaccadeDurationsMillisecondsList = ...
        obj.microSaccadeDurationsMillisecondsDist(1:maxMicroSaccadeNum);
end

% Zero the microSaccadeResidualPath
obj.microSaccadeResidualPath = zeros(2,1);

end
