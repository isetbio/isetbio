classdef fixationalEM < handle
%FIXATIONALEM Class - generates fixational eye movements and saccades
%
% Description:
%   The drift EM component is generated using the model described in: 
%   Mergenthaler and Engbert (2007) Modeling the control of fixational eye movements with
%   neurophysiological delays", Physical Review Letters, 98 (13)
%
% Syntax:
%   em = fixationalEM();
%
% Inputs:
%   None
% Returns:
%   em = Fixational EM object
%
% Optional key/value pairs
%    
% See also
%   emOBJ.compute(durationSeconds, sampleTimeSeconds);
%
% Public properties (user-settable)
properties    
    controlNoiseMean;
    controlNoiseStd;       
    controlGamma;              
    feedbackGain;             
    feedbackSteepness;          
    feedbackXposDelaySeconds;
    feedbackYposDelaySeconds;
    positionNoiseMean;
    positionNoiseStd; 
    
    stabilizationSeconds;                           % how long to warm up the drift model for

    % Microsaccade-related properties
    microSaccadeType;                               % choose between {'none', 'heatmap/fixation based', 'stats based'}
    microSaccadeMeanIntervalSeconds;                % mean of the Gamma distribution of microsaccade intervals
    microSaccadeMeanAmplitudeArcMin;                % mean of the Gamma distribution of microsaccade amplitudes
    microSaccadeMeanSpeedDegsPerSecond;             % mean of the Normal distribution of microsaccade speeds
    microSaccadeStDevSpeedDegsPerSecond;            % std of the Normal distribution of microsaccade speeds
    microSaccadeTargetJitterArcMin;                 % jitter in microsaccade target
    microSaccadeDirectionJitterDegs;                % jitter in microsaccade direction (only used in stats-only mode)
    microSaccadeMinDurationMilliSeconds;            % minumum duration of a microsaccade

    % Fixation map properties
    fixationMapSpaceConstantArcMin;                 % sigma of Gaussian specifying likelihood of a saccade position (from 0,0)
   
    % Heat map properties
    heatMapWeight;                                  % How much to weight the heat map vs. the fixation map
    heatMapWidthArcMin;                             % Width of heat map
    heatMapSpatialSampleArcMin;                     % Spatial resolution of heat map
    heatMapTemporalSampleSeconds;                   % Temporal resolution of heat map
    heatMapKernelSpaceConstantArcMin;               % Space constant for heat map's spatial kernel
    heatMapKernelTimeConstantSeconds;               % Time constant for heat map's temporal kernel
    
    randomSeed;                                     % Random seed to use
    
    % user-feedback
    displayComputeProgress;                         % Set to true, to display a progress bar
    beVerbose;                                      % Set to false, no printing
end

% Read-only properties
properties (SetAccess = private)   
    timeAxis = [];                                  % time axis of the em path
    velocityArcMin = [];                            % nTrials x emPath length     matrix of velocity in units of arc minutes/second
    emPosArcMin = [];                               % nTrials x emPath length x 2 matrix of eye positions in units of arc minutes - no truncation
    emPosMicrons  = [];                             % nTrials x emPath length x 2 matrix of eye positions in units of microns - truncated to cone mosaic pattern size
    emPos = [];                                     % nTrials x emPath length x 2 matric of eye positions in units of the cone mosaic pattern size
end

properties (GetAccess = protected)
    % The properties  contain data that are usefull for visualizing (e.g., demo app)
    % Their values can be obtained via the subclass method
    % @visualizableFixationalEM.getValueOfProperty(propertyName)
    
    emPosTimeSeriesArcMin = [];                     % emPath - last trial
    velocityArcMinPerSecTimeSeries = [];            % velocity - last trial
    
    microSaccadeIntervalsMilliSecondsDist = [];     % distribution of microsaccade intervals - useful for visualization
    microSaccadeIntervalsMilliSecondsList = [];     % list with possible microsaccade intervals
    microSaccadeIntervalsMilliseconds = []          % list with generated microsaccade intervals
    
    microSaccadeAmplitudesArcMinDist = [];          % distribution of microsaccade amplitudes - useful for visualization
    microSaccadeAmplitudesArcMinList = [];          % list with possible microsaccade amplitudes
    microSaccadeAmplitudesArcMin = [];              % list with generated microsaccade amplitudes

    microSaccadeSpeedDegsPerSecondDist = [];        % distribution of microsaccade speeds - useful for visualization
    microSaccadeSpeedDegsPerSecondList = [];        % list with possible microsaccade speeds
    microSaccadeSpeedsDegsPerSecond = [];           % list with generated microsaccade speeds
    
    microSaccadeDurationsMillisecondsDist = [];     % distribution of microsaccde durations (derived from speed and amplitude distribution)
    microSaccadeDurationsMillisecondsList = [];     % list with generated microsaccade durations
    
    microSaccadeOnsetStepIndices = [];              % list with generated microsaccade onset time steps
    microSaccadeTargetLikelihoodMapsList = [];      % list with generated microsaccade target likelihood maps
    microSaccadePositions = [];
     
    heatMapTimeSeries = [];                         % The heat map generated by the emPath
    heatMapXYsupport = [];                          % Spatial support for the heatmap
end


% Properties invisible to the user
properties (Access = private)
    stabilizationStepsNum = [];      
    tStepsNum = [];                                 % total number of time steps for the entire duration of the drift
    feedbackXposDelayTimeSteps;                     % feedback delay for the x-pos in time steps
    feedbackYposDelayTimeSteps;                     % feedback delay for the y-pos in time steps
    
    heatMapTimeSeriesIntermediate = []; 
    emAccumMap = [];
    heatMapUpdateIntervalStepsNum = [];
    lastMicroSaccadeTimeStep = [];
    emPosTimeSeries= [];                            % time series of emPositions (preallocated in initDriftState)
    velocityTimeSeries = [];
    
    % Microsaccade - related
    microSaccadeIntervalGammaShapeParameter;        % Shape parameter of the Gamma distribution of microsaccade intervals
    microSaccadeAmplitudeGammaShapeParameter;       % Shape parameter of the Gamma distribution of microsaccade amplitudes
    microSaccadeDurationGammaShapeParameter;        % Shape parameter of the Gamma distribution of microsaccade durations
    microSaccadeTargetLikelihoodSpatialMap = [];    % The likelihood map for saccades (Gaussian centered at (0,0)) with sigma : obj.fixationMapSpaceConstantArcMin
    microSaccadeResidualPath = [];                  % multi-step microsaccade path
    validMicrosaccadeAmplitudeIndices = [];         % list with amplitudes that have been already used (sample with no replacement)
    
    % Drift model internal variables time series
    positionalNoiseTimeSeries = [];                 % time series of positional noise  (precomputed in initDriftState) - last trial
    controlNoiseTimeSeries = [];                    % time series of control signal noise (precomputed in initDriftState)
    controlSignalTimeSeries = [];                   % time series of the control signal  (preallocated in initDriftState)
    feedbackSignalTimeSeries = [];                  % time series of the feedback signal (preallocated in initDriftState)

    currentPositionHeatLevelTimeSeries = [];        % The heat level at the current position - could be used to trigger a saccade, but currently not used
    heatMapSpatialKernel = [];
    heatMapSpatialKernelSupport = [];
    heatMapTemporalKernel = [];
    heatMapTemporalKernelSupport = [];
end

properties (Constant)
    timeStepDurationSeconds = 1/1000;
    velocityMeasurementIntervalSeconds = 41/1000;  % window for measuring velocity to 41 milliseconds, as in Cherici & Rucci 2012 - Precision of sustained fixation in trained and untrained observers)
    scalarToArcMin = 0.5;                          % a 0.5 scaling factor gives drift speed in the middle of the range reported byCherici et al 2012
end

methods
    % Constructor
    function obj = fixationalEM(varargin)
        setDefaultParams(obj);
    end
    
    % Method
	setDriftParamsFromMergenthalerAndEngbert2007Paper(obj);

    % Method
    setMicroSaccadeStats(obj);

    % Method
    setDefaultParams(obj);
    
    % Method
	computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, varargin);
    
    % Method
    compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials, computeVelocity, varargin);

    % Method to compute an emPath with duration emDurationSeconds
    % sampled with a temporal resolution of sampleDurationSeconds
    computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds);
end % Public methods

methods (Static)
    % Method
    [fixationMap, fixationMapSupportX, fixationMapSupportY, fixationMapXSlice, fixationMapYSlice] = ...
        computeFixationMap(timeAxis, emPaths, emPosRange, emPosDelta, varargin);

    % Method
    [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        performDisplacementAnalysis(emPos, timeAxisSeconds, varargin);

    % Method
    [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        performD2DisplacementAnalysis(emPos, timeAxisSeconds);
    
    % Method to generate a video of fixationalEMs together with a cone mosaic
    generateEMandMosaicComboVideo(fixationalEMObj, coneMosaicHexObj, varargin);
end % Static methods

methods (Access = private)

    % Method
    initOutputs(obj);
    
    % Method to initiate the state of the object
    initState(obj, emDurationSeconds);

    % Method
    checkForMicroSaccadeEpoch(obj,tStep);

    % Method
    updateMicroSaccadeResidualPath(obj);

    % Method
    targetPosArcMin = posOfMaxLikelihood(obj,targetLikelihoodSpatialMap);

    % Method
    currentPositionHeatLevel = updateHeatMap(obj, tStep);

    % Method
    velocityTimeSeries = computeVelocity(obj, theEMpath);

    % Method
    trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds);

    % Method
    outMatrix = smartInterpolation(obj,inputTimeAxis, inputMatrix, outputTimeAxis);

    % Method to compute the spatial support for the heatMap
    heatMapXYsupport = spatialSupportForHeatMap(obj);

    % Method to compute the spatial kernel for the heat map
    [kernel, kernelSupport] = spatialKernelForHeatMap(obj);

    % Method to compute the temporal kernel for the heat map
    [kernel, kernelSupport] = temporalKernelForHeatMap(obj, gridNodesNum);

    % Method to compute obj.microSaccadeTargetLikelihoodSpatialMap
    microSaccadeTargetLikelihoodSpatialMap = computeMicrosaccadeTargetLikelihoodSpatialMap(obj, heatMapXYsupportArcMin);
end % private methods

end % class