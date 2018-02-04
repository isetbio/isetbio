classdef fixationalEM < handle
%FIXATIONALEM Class for generating fixational eye movements using realistic drift
% and saccade models
%
% Syntax:
%   emOBJ = fixationalEM();
%   emOBJ.compute(durationSeconds, sampleTimeSeconds);
%
% Description:
%   The drift EM component is generated using the model described in: 
%   Mergenthaler and Engbert (2007) Modeling the control of fixational eye movements with
%   neurophysiological delays", Physical Review Letters, 98 (13)

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
    velocityMeasurementIntervalSeconds;             % interval over which to measure velocity

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
    timeAxis = [];              % time axis of the em path
    emPosArcMin = [];           % nTrials x emPath length x 2 matrix of eye positions in units of arc minutes - no truncation
    emPosMicrons  = [];         % nTrials x emPath length x 2 matrix of eye positions in units of microns - truncated to cone mosaic pattern size
    emPos = [];                 % nTrials x emPath length x 2 matric of eye positions in units of the cone mosaic pattern size
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
    tStepsNum = [];                  % total number of time steps for the entire duration of the drift
    feedbackXposDelayTimeSteps;      % feedback delay for the x-pos in time steps
    feedbackYposDelayTimeSteps;      % feedback delay for the y-pos in time steps
    
    heatMapTimeSeriesIntermediate = []; 
    emAccumMap = [];
    heatMapUpdateIntervalStepsNum = [];
    lastMicroSaccadeTimeStep = [];
    emPosTimeSeries= [];             % time series of emPositions (preallocated in initDriftState)
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
    scalarToArcMin = 1.0;            % factor to scale positions to arc min
end

methods
% Constructor
function obj = fixationalEM(varargin)
    setDefaultParams(obj);
end

function setDriftParamsFromMergenthalerAndEngbert2007Paper(obj)
    % Control gain (gamma in Mergenthaler&Engbert, 2007)
    obj.controlGamma = 0.25;
    
    % Control noise, ksi, with zero MEAN and sigma STD
    ksi = struct('mean', 0, 'sigma', 0.075);
    obj.controlNoiseMean = ksi.mean;
    obj.controlNoiseStd = ksi.sigma; 
    
    % Position noise, eta with zero MEAN and rho STD
    eta = struct('mean', 0, 'rho', 0.35/2);
    obj.positionNoiseMean = eta.mean;
    obj.positionNoiseStd = eta.rho;    
    
    %Feedback gain (lambda in Mergenthaler&Engbert, 2007)
    obj.feedbackGain = 0.15;             
     
    % Feedback steepness (epsilon in Mergenthaler&Engbert, 2007)
    obj.feedbackSteepness = 1.1;          
    
    % Feedback delay for x- and y-pos 
    obj.feedbackXposDelaySeconds = 70/1000; % tau_hor in Mergenthaler&Engbert, 2007
    obj.feedbackYposDelaySeconds = 40/1000; % tau_vert in Mergenthaler&Engbert, 2007
end

function setMicroSaccadeStats(obj)
    % Set the inter-saccde interval distribution
    obj.microSaccadeMeanIntervalSeconds = 0.45;
    obj.microSaccadeIntervalGammaShapeParameter = 5;
    
    % Set the microsaccade amplitude distribution 
    % (consistent with Fig 1 in Martinez 2009)
    obj.microSaccadeMeanAmplitudeArcMin = 8;
    obj.microSaccadeAmplitudeGammaShapeParameter = 4;
    
    % 30 deg/second (See Martinez 2008, Table in Fig 3D)
    obj.microSaccadeMeanSpeedDegsPerSecond = 39;
    obj.microSaccadeStDevSpeedDegsPerSecond = 2;  
    
    % Minimum duration of a microsaccade - Arbitrary value
    obj.microSaccadeMinDurationMilliSeconds = 2;
    
    % Sigma of microsaccade target jitter - Arbitrary value
    obj.microSaccadeTargetJitterArcMin = 0.3;
    
    % Sigma of corrective microsaccade direction jitter - Arbitrary value
    obj.microSaccadeDirectionJitterDegs = 15;
end

function setDefaultParams(obj)
    
    % Set the drift params from the Mergenthaler&Engbert, 2007 paper
    obj.setDriftParamsFromMergenthalerAndEngbert2007Paper();
    
    % Set microsaccade stats consistent with data in Martinez, 2009 
    obj.setMicroSaccadeStats();
    
    % Set the window for measuring velocity to 41 milliseconds
    % (as in Cherici, Rucci 2012 - Precision of sustained fixation in trained and untrained observers)
    obj.velocityMeasurementIntervalSeconds = 41/1000;

    % How long to 'warm' up the drift model for
    obj.stabilizationSeconds = 2.0;
    
    % Microsaccade generation - related params
    obj.microSaccadeType = 'heatmap/fixation based';
    obj.heatMapWeight = 0.5;
    
    % Fixation map properties
    obj.fixationMapSpaceConstantArcMin = 15;

    % Heat map properties
    obj.heatMapWidthArcMin = 50;
    obj.heatMapSpatialSampleArcMin = 2;
    obj.heatMapTemporalSampleSeconds = 10/1000;
    obj.heatMapKernelSpaceConstantArcMin = 2.0;
    obj.heatMapKernelTimeConstantSeconds = 0.5;
    
    % Do not set the random seed
    obj.randomSeed = [];
    
    % Compute progress display
    obj.displayComputeProgress = false;
    obj.beVerbose = false;
end

function computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, varargin)
    p = inputParser;
    p.addRequired('coneMosaic', @(x)(isa(x, 'coneMosaic')));
    p.addRequired('eyeMovementsPerTrial', @isscalar);
    p.addParameter('nTrials',1,@isscalar);
    p.addParameter('rSeed', [], @isscalar);
    p.parse(coneMosaic, eyeMovementsPerTrial, varargin{:});
    
    % Set optional parameters based on input
    obj.randomSeed = p.Results.rSeed;
    nTrials = p.Results.nTrials;
    
    % Extract sampleDuration from coneMosaic's integration time
    sampleDurationSeconds = coneMosaic.integrationTime;
    
    % Compute emDurationSeconds
    emDurationSeconds = eyeMovementsPerTrial * sampleDurationSeconds;
    
    if (isempty(obj.randomSeed))
        rng('shuffle');
    else
        rng(obj.randomSeed);
    end
    
    % Generate fixational eye movements (obj.emPathsArcMin)
    compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials);
    
    % Subsample to cone positions
    conePatternSampleMicrons = coneMosaic.patternSampleSize(1)*1e6;
    emPosMicronsNotRounded = (obj.emPosArcMin/60)*coneMosaic.micronsPerDegree;
    % Return the paths in units of cone mosaic pattern size
    obj.emPos = sign(emPosMicronsNotRounded) .* round(abs(emPosMicronsNotRounded)/conePatternSampleMicrons);
    % Also return the path in units of microns
    obj.emPosMicrons = obj.emPos * conePatternSampleMicrons;
end

function compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials)
    
    if (isempty(obj.randomSeed))
        rng('shuffle');
    else
        rng(obj.randomSeed);
    end
    
    for iTrial = 1:nTrials
        if (obj.beVerbose)
            fprintf('Computing emModel for %2.2f seconds; sampleT = %2.4f sec, samples: %d\n', emDurationSeconds, sampleDurationSeconds, obj.tStepsNum);
        end
        computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds);
        if (iTrial == 1)
            obj.emPosArcMin = zeros(nTrials, length(obj.emPosTimeSeriesArcMin), 2);
        end
        obj.emPosArcMin(iTrial,:,:) = reshape(obj.emPosTimeSeriesArcMin', [1 length(obj.emPosTimeSeriesArcMin) 2]);
    end
end

% Method to compute an emPath with duration emDurationSeconds
% sampled with a temporal resolution of sampleDurationSeconds
function computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds)
    
    % Init the state
    initState(obj, emDurationSeconds);
    
    feedbackDelaySteps(1) = obj.feedbackXposDelayTimeSteps;
    feedbackDelaySteps(2) = obj.feedbackYposDelayTimeSteps;

    % form time axis
    obj.timeAxis = ((1:obj.tStepsNum)-1) * obj.timeStepDurationSeconds;

    % Init the time step when the last saccade occurred
    obj.lastMicroSaccadeTimeStep = obj.stabilizationStepsNum;
         
    if (obj.displayComputeProgress)
        hWaitBar = waitbar(0, 'Computing emPath ...');
    end 
    
    for tStep = 1:obj.tStepsNum
        
        if (mod(tStep-1, 100) == 0) && (obj.displayComputeProgress)
            progress = tStep/obj.tStepsNum;
            waitbar(progress, hWaitBar, 'Computing emPath ...');
        end
        for eyeIndex = 1:2
            % compute feedback signals
            if (tStep > feedbackDelaySteps)
                delayedControl = obj.controlSignalTimeSeries(eyeIndex, tStep-feedbackDelaySteps(eyeIndex));
                obj.feedbackSignalTimeSeries(eyeIndex,tStep) = ...
                    obj.feedbackGain * tanh(obj.feedbackSteepness*delayedControl);
            end

            % compute next control signals
            obj.controlSignalTimeSeries(eyeIndex,tStep+1) = ...
                (1-obj.controlGamma)*obj.controlSignalTimeSeries(eyeIndex,tStep) - ...
                obj.feedbackSignalTimeSeries(eyeIndex,tStep) + obj.controlNoiseTimeSeries(eyeIndex,tStep);

            % generate next EMpoint
            obj.emPosTimeSeries(eyeIndex,tStep+1) = ...
                obj.emPosTimeSeries(eyeIndex,tStep) + ...
                obj.controlSignalTimeSeries(eyeIndex,tStep+1) + ...
                obj.positionalNoiseTimeSeries(eyeIndex, tStep) + ...
                obj.microSaccadeResidualPath(eyeIndex,1);
        end % eyeIndex

        % Update the microSaccadeResidualPath
        obj.updateMicroSaccadeResidualPath();
            
        % Check for saccade generation only after the drift has stabilized
        if (tStep>=obj.stabilizationStepsNum)
            if strcmp(obj.microSaccadeType ,'heatmap/fixation based')
                % Update heat map and retrieve heat level at the current emPos. 
                % Maybe we could also trigger saccades when the
                % currentPositionHeatLevel exceeds a threshold.
                currentPositionHeatLevel = obj.updateHeatMap(tStep);
            end
            if (~strcmp(obj.microSaccadeType,'none'))
                % Check whether a microsaccade will be generated
                obj.checkForMicroSaccadeEpoch(tStep);
            end
        end % if (tStep>=obj.stabilizationStepsNum)
    end % tStep

    % Compute velocity
    obj.velocityTimeSeries = obj.computeVelocity(obj.timeAxis, obj.emPosTimeSeries', obj.velocityMeasurementIntervalSeconds);
    
    % Trim and recenter and resample the post-stabilization time series
    trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds);

    if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
        % Scale heat map to 1 (this affects only visualization at this point)
        obj.heatMapTimeSeries = obj.heatMapTimeSeries / max(obj.heatMapTimeSeries(:));
    end
    
    % Make the arc min sequence
    obj.emPosTimeSeriesArcMin = obj.emPosTimeSeries*obj.scalarToArcMin;   
    obj.velocityArcMinPerSecTimeSeries = obj.velocityTimeSeries*obj.scalarToArcMin;
    
    if (obj.displayComputeProgress) && (obj.displayComputeProgress)
        close(hWaitBar);
    end
end % compute 

end % Public methods

methods (Static)
    
function  velocityTimeSeries = computeVelocity(timeAxis, emPosTimeSeries, velocityMeasurementIntervalSeconds)
    emPosTimeSeries = emPosTimeSeries';
    timeStepDurationSeconds = timeAxis(2)-timeAxis(1);
    method = 2;
    if (method == 1)
        velocityMeasurementIntervalHalfSteps = floor(velocityMeasurementIntervalSeconds/timeStepDurationSeconds/2);
        sigmaSteps = velocityMeasurementIntervalHalfSteps/2.5;
        tSteps = (-velocityMeasurementIntervalHalfSteps:velocityMeasurementIntervalHalfSteps);
        weightKernel = exp(-0.5*(tSteps/sigmaSteps).^2);
        weightKernel = weightKernel / sum(weightKernel(:));
        smoothedPos = conv2(emPosTimeSeries, weightKernel, 'same');
        diffSmoothedPos = diff(smoothedPos, 1, 2);
        velocities = conv2(diffSmoothedPos, weightKernel, 'same')/timeStepDurationSeconds;
        velocityTimeSeries = sqrt(squeeze(sum(velocities.^2, 1)));
    else
        % make smoothing kernel
        velocityIntervalHalfSamples = max([1 floor(velocityMeasurementIntervalSeconds/timeStepDurationSeconds/2)]);
        stepsToAverage = 2*velocityIntervalHalfSamples;
        dd = -stepsToAverage:stepsToAverage;
        sigmaD = stepsToAverage/3;
        weightKernel = exp(-0.5*(dd'/sigmaD).^2);
        weightKernel = weightKernel/sum(weightKernel);
        weightKernel = [weightKernel(:) weightKernel(:)]';
    
        meanDeltas = nan(1,length(emPosTimeSeries));
        firstTsample = 1+stepsToAverage+velocityIntervalHalfSamples;
        lastTsample = length(emPosTimeSeries)-stepsToAverage-velocityIntervalHalfSamples;
        
        for t = firstTsample:lastTsample
            deltas = abs(emPosTimeSeries(:,t+dd+velocityIntervalHalfSamples)-emPosTimeSeries(:,t+dd-velocityIntervalHalfSamples));
            deltas = sum(deltas .* weightKernel, 2);
            meanDeltas(t) = sqrt(sum(deltas.^2));
        end
        velocityTimeSeries = meanDeltas/((2*velocityIntervalHalfSamples+1)*timeStepDurationSeconds);
    end
end

function [emLikelihoodMap, emLikelihoodMapSupportX, emLikelihoodMapSupportY] = computeLikelihoodMap(emPaths, emPosRange, emPosDelta)
    xBins = emPosRange(1):emPosDelta:emPosRange(end)+emPosDelta;
    yBins = emPosRange(1):emPosDelta:emPosRange(end)+emPosDelta;
    nTrials = size(emPaths,1);
    for trialNo = 1:nTrials
        emPath = squeeze(emPaths(trialNo,:,:));
        [N, xEdges, yEdges] = histcounts2(squeeze(emPath(:,1)), squeeze(emPath(:,2)), xBins, yBins);
        if (trialNo == 1)
            emLikelihoodMap = N;
        else
            emLikelihoodMap = emLikelihoodMap + N;
        end
    end
    % Map across all trials
    emLikelihoodMap = emLikelihoodMap/max(emLikelihoodMap(:));
    emLikelihoodMapSupportX = xEdges(1:end-1)+emPosDelta/2;
    emLikelihoodMapSupportY = yEdges(1:end-1)+emPosDelta/2;
end

function [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        performDisplacementAnalysis(emPos, timeAxisSeconds, varargin)
%
    p = inputParser;
    addParameter(p, 'mode', 'D1', @(x) (ischar(x) && ismember(x, {'D1', 'D2'})));
    parse(p, varargin{:});
    
    timeLagsMilliseconds = [1 2:2:100 105:5:1000];
    dtMilliseconds = (timeAxisSeconds(2)-timeAxisSeconds(1))*1000;

    % Compute emPos from randomized intervals
    sdiff = diff(emPos);
    scrambledIndices = randperm(length(sdiff));
    sdiff = sdiff(scrambledIndices);
    emPosScrambled = [emPos(1) emPos(1)+cumsum(sdiff)];

    % Initialize displacement vectors
    meanD = zeros(1,numel(timeLagsMilliseconds));
    maxD = zeros(1,numel(timeLagsMilliseconds));
    meanDscrambled = zeros(1,numel(timeLagsMilliseconds));

    tStepsNum = length(timeAxisSeconds);
    for delayIndex = 1:numel(timeLagsMilliseconds)
        m = round(timeLagsMilliseconds(delayIndex)/dtMilliseconds);
        sumDeltaRsquared = zeros(1, tStepsNum-m);
        sumDeltaRsquaredScrambled = zeros(1, tStepsNum-m);
        for i = 1:tStepsNum-m 
            sumDeltaRsquared(i) = abs(emPos(i+m)-emPos(i));
            sumDeltaRsquaredScrambled(i) = abs(emPosScrambled(i+m)-emPosScrambled(i));
        end
        if (strcmp(p.Results.mode, 'D2'))
            sumDeltaRsquared = sumDeltaRsquared.^2;
            sumDeltaRsquaredScrambled = sumDeltaRsquaredScrambled.^2;
        end
        if (numel(sumDeltaRsquared) > 0)
            maxD(1,delayIndex) = max(sumDeltaRsquared);
            meanD(1,delayIndex) = mean(sumDeltaRsquared);
            meanDscrambled(1,delayIndex) = mean(sumDeltaRsquaredScrambled);
        end
    end
    k = meanDscrambled(1)/timeLagsMilliseconds(1);
    H = 0.5;
    DrandomWalk = k*timeLagsMilliseconds.^H;
end

function [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = performD2DisplacementAnalysis(emPos, timeAxisSeconds)
    timeLagsMilliseconds = [1 2:2:100 105:5:1000];
    dtMilliseconds = (timeAxisSeconds(2)-timeAxisSeconds(1))*1000;

    % Compute emPos from randomized intervals
    sdiff = diff(emPos);
    scrambledIndices = randperm(length(sdiff));
    sdiff = sdiff(scrambledIndices);
    emPosScrambled = [emPos(1) emPos(1)+cumsum(sdiff)];

    % Initialize displacement vectors
    meanD = nan(1,numel(timeLagsMilliseconds));
    maxD = nan(1,numel(timeLagsMilliseconds));
    meanDscrambled = nan(1,numel(timeLagsMilliseconds));

    tStepsNum = length(timeAxisSeconds);
    for delayIndex = 1:numel(timeLagsMilliseconds)
        m = round(timeLagsMilliseconds(delayIndex)/dtMilliseconds);
        sumDeltaRsquared = zeros(1, tStepsNum-m);
        sumDeltaRsquaredScrambled = zeros(1, tStepsNum-m);
        for i = 1:tStepsNum-m 
            sumDeltaRsquared(i) = (norm(emPos(i+m)-emPos(i)))^2;
            sumDeltaRsquaredScrambled(i) = (norm(emPosScrambled(i+m)-emPosScrambled(i)))^2;
        end
        if (numel(sumDeltaRsquared) > 0)
            maxD(delayIndex) = max(sumDeltaRsquared);
            meanD(delayIndex) = mean(sumDeltaRsquared);
            meanDscrambled(delayIndex) = mean(sumDeltaRsquaredScrambled);
        end
    end
    k = meanDscrambled(1)/timeLagsMilliseconds(1);
    H = 0.5;
    DrandomWalk = k*timeLagsMilliseconds.^H;
end

end % Static methods

methods (Access = private)
% Init the state of the object
function initState(obj, emDurationSeconds)     
    % Compute stabilization steps
    obj.stabilizationStepsNum = round(obj.stabilizationSeconds / obj.timeStepDurationSeconds);

    % Compute total number of time steps 
    obj.tStepsNum = round(emDurationSeconds / obj.timeStepDurationSeconds) + obj.stabilizationStepsNum + 1;

    % Pre-allocate memory for time axis                
    obj.timeAxis = zeros(1, obj.tStepsNum);

    % Drift - related params
    % Compute feedback delay in time steps for X/Y position
    obj.feedbackXposDelayTimeSteps = round(obj.feedbackXposDelaySeconds/obj.timeStepDurationSeconds);
    obj.feedbackYposDelayTimeSteps = round(obj.feedbackYposDelaySeconds/obj.timeStepDurationSeconds);

    % Compute positional noise for the entire duration of the drift
    obj.positionalNoiseTimeSeries = normrnd(...
        obj.positionNoiseMean * ones(2, obj.tStepsNum), ...
        obj.positionNoiseStd  * ones(2, obj.tStepsNum));

    % Compute control signal noise for the entire duration of the drift
    obj.controlNoiseTimeSeries = normrnd(...
        obj.controlNoiseMean * ones(2, obj.tStepsNum), ...
        obj.controlNoiseStd  * ones(2, obj.tStepsNum));

    % Pre-allocate memory for the feedback and the control components of the drift mechanism
    obj.feedbackSignalTimeSeries = zeros(2, obj.tStepsNum);
    obj.controlSignalTimeSeries = zeros(2, obj.tStepsNum);

    % Pre-allocate memory for the computed emPositions - assuming 1 millisecond time sample
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
    obj.microSaccadeTargetLikelihoodSpatialMap = obj.computeMicrosaccadeTargetLikelihoodSpatialMap(obj.heatMapXYsupport);
        
    tBins = numel(obj.stabilizationStepsNum:obj.tStepsNum);
    obj.heatMapTimeSeries = zeros(tBins, numel(obj.heatMapXYsupport), numel(obj.heatMapXYsupport));
    obj.currentPositionHeatLevelTimeSeries = zeros(1, tBins);

    obj.heatMapUpdateIntervalStepsNum = round(obj.heatMapTemporalSampleSeconds/obj.timeStepDurationSeconds);
    tBinsSubSampled = 1:obj.heatMapUpdateIntervalStepsNum:numel(obj.stabilizationStepsNum);
    obj.heatMapTimeSeriesIntermediate = zeros(numel(tBinsSubSampled), numel(obj.heatMapXYsupport), numel(obj.heatMapXYsupport));

    % Compute spatial kernel for heat map
    [obj.heatMapSpatialKernel, obj.heatMapSpatialKernelSupport] = obj.spatialKernelForHeatMap();

    % Compute temporal kernel for heat map
    [obj.heatMapTemporalKernel, obj.heatMapTemporalKernelSupport] = ...
        obj.temporalKernelForHeatMap(length(obj.heatMapXYsupport));
    
    % Generate Gamma distibuted microsaccade intervals
    maxMicroSaccadeNum = obj.tStepsNum - obj.stabilizationStepsNum;
    shape = obj.microSaccadeIntervalGammaShapeParameter; 
    scale = obj.microSaccadeMeanIntervalSeconds*1000/shape;
    obj.microSaccadeIntervalsMilliSecondsDist = round(random('Gamma', shape, scale, [1 max([10000 obj.tStepsNum])]));
    obj.microSaccadeIntervalsMilliSecondsList = obj.microSaccadeIntervalsMilliSecondsDist(1:maxMicroSaccadeNum);
    
    % Generate microsaccade speed distribution
    obj.microSaccadeSpeedDegsPerSecondDist = random('Normal', obj.microSaccadeMeanSpeedDegsPerSecond, obj.microSaccadeStDevSpeedDegsPerSecond, [1 max([10000 obj.tStepsNum])]);
    obj.microSaccadeSpeedDegsPerSecondList = obj.microSaccadeSpeedDegsPerSecondDist(1:maxMicroSaccadeNum);
        
    if (strcmp(obj.microSaccadeType, 'stats based'))
        % Generate Gamma distributed microsaccade amplitudes
        shape = obj.microSaccadeAmplitudeGammaShapeParameter; 
        scale = obj.microSaccadeMeanAmplitudeArcMin/shape;
        obj.microSaccadeAmplitudesArcMinDist = random('Gamma', shape, scale, [1 max([10000 obj.tStepsNum])]);
        obj.microSaccadeAmplitudesArcMinList = obj.microSaccadeAmplitudesArcMinDist(1:maxMicroSaccadeNum);
        obj.validMicrosaccadeAmplitudeIndices = 1:numel(obj.microSaccadeAmplitudesArcMinList);
        
        % Generate microsaccade durations distribution from the speed and the amplitude distributions
        obj.microSaccadeDurationsMillisecondsDist = round(obj.microSaccadeAmplitudesArcMinDist./(obj.microSaccadeSpeedDegsPerSecondDist/1000*60));
        obj.microSaccadeDurationsMillisecondsDist(obj.microSaccadeDurationsMillisecondsDist<obj.microSaccadeMinDurationMilliSeconds) = obj.microSaccadeMinDurationMilliSeconds;
        obj.microSaccadeDurationsMillisecondsList = obj.microSaccadeDurationsMillisecondsDist(1:maxMicroSaccadeNum);
    end
    
    % Zero the microSaccadeResidualPath
    obj.microSaccadeResidualPath = zeros(2,1);
end % initState(obj)


function checkForMicroSaccadeEpoch(obj,tStep)

    % Compute probablility of generating a saccade base on interval since last saccade
    intervalSinceLastMicroSaccadeSeconds = (tStep-obj.lastMicroSaccadeTimeStep)*obj.timeStepDurationSeconds;
    
    if (round(intervalSinceLastMicroSaccadeSeconds*1000) ~= obj.microSaccadeIntervalsMilliSecondsList(numel(obj.microSaccadeOnsetStepIndices)+1))
        % Not generating a microsaccade this time
        return;
    end
   
    % Ok, we are generating a saccade
    % Retrieve current position (aligned with respect to eye position at the end of the stabilization time)
    currentPosCentered = obj.emPosTimeSeries(:,tStep+1)-obj.emPosTimeSeries(:,obj.stabilizationStepsNum);
    
    % Generate a positional jitter vector
    shape = 5;
    scale = obj.microSaccadeTargetJitterArcMin/shape;
    positionJitter = random('Gamma', shape, scale, [1 1]); % randn(2,1).*obj.microSaccadeTargetJitterArcMin
    positionJitterNorm = norm(positionJitter);
    
    % Determine microsaccade target. 
    if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
        % Compute normalized currentHeatMap
        postStabilizationStep = 1+tStep-obj.stabilizationStepsNum;
        currentHeatMap = squeeze(obj.heatMapTimeSeries(postStabilizationStep,:,:));
        currentHeatMap = currentHeatMap/max(currentHeatMap(:));
        % Compute combined saccade target likelihood spatial map
        targetLikelihoodSpatialMap = obj.heatMapWeight*(1-currentHeatMap) + (1-obj.heatMapWeight)*obj.microSaccadeTargetLikelihoodSpatialMap;
        % Find point at which the targetLikelihoodSpatialMap is maximal
        targetPosArcMin = obj.posOfMaxLikelihood(targetLikelihoodSpatialMap);
        targetPos = targetPosArcMin/obj.scalarToArcMin;
        % Compute the displacement vector
        microSaccadeJumpArcMin = (targetPos'-currentPosCentered)*obj.scalarToArcMin + positionJitter;
        % Compute the duration based on the displacement vector and a speed from the obj.microSaccadeSpeedDegsPerSecondList
        saccadeIndex = numel(obj.microSaccadeOnsetStepIndices)+1;
        speedArcMinPerMillisecond = obj.microSaccadeSpeedDegsPerSecondList(saccadeIndex) / 1000*60;
        durationMilliSeconds = round(norm(microSaccadeJumpArcMin) / speedArcMinPerMillisecond); 
        durationMilliSeconds = max([obj.microSaccadeMinDurationMilliSeconds durationMilliSeconds]);
        
    elseif (strcmp(obj.microSaccadeType, 'stats based'))
        % Decide whether we will do a corrective saccade based on how far away we are 
        correctionDecisionFactor = 1-exp(-(norm(currentPosCentered*obj.scalarToArcMin)/(0.5*obj.fixationMapSpaceConstantArcMin))^2);
        doCorrectiveSaccade = correctionDecisionFactor > rand(1,1);
        if (doCorrectiveSaccade)
            % Corrective saccade
            % direction toward the fixation point
            theta = atan2(currentPosCentered(2),currentPosCentered(1)) + pi + randn(1,1)*obj.microSaccadeDirectionJitterDegs/180*pi;
            
            % Distance to fixation
            distanceToFixationPointArcMin = norm(squeeze(currentPosCentered)) * obj.scalarToArcMin;
            % Add some jitter
            distanceToFixationPointArcMin = distanceToFixationPointArcMin + positionJitterNorm;
    
            % Find the index of closest amplitude from the obj.microSaccadeAmplitudesArcMinList
            [~,idx] = min(abs(distanceToFixationPointArcMin-obj.microSaccadeAmplitudesArcMinList(obj.validMicrosaccadeAmplitudeIndices)));
            idx = obj.validMicrosaccadeAmplitudeIndices(idx);
        else
            % Random saccade
            % random direction, random amplitude, but one that does not take us further away than some threshold distance
            theta = rand(1,1)*2*pi;
            
            % Find an amplitude index that is within the obj.fixationMapSpaceConstantArcMin region
            possibleJumps = (obj.microSaccadeAmplitudesArcMinList(obj.validMicrosaccadeAmplitudeIndices))' * [cos(theta) sin(theta)];
            possibleEndPoints = bsxfun(@plus, possibleJumps, (currentPosCentered')*obj.scalarToArcMin);
            possibleDistancesToFixationPoint = sqrt(sum(possibleEndPoints.^2, 2));
            idx = find(possibleDistancesToFixationPoint < obj.fixationMapSpaceConstantArcMin);
            if (isempty(idx))
                [~,idx] = min(possibleDistancesToFixationPoint);
            else
                randomIdx = randperm(numel(idx));
                idx = idx(randomIdx(1));
            end
            idx = obj.validMicrosaccadeAmplitudeIndices(idx);
        end
        
        % Retrieve the amplitude at the computed index
        amplitudeArcMin = obj.microSaccadeAmplitudesArcMinList(idx);
    
        % Retrieve the correspondng duration
        durationMilliSeconds = obj.microSaccadeDurationsMillisecondsList(idx);
        
        % Remove this amplitude from the list of valid indices (sampling without replacement)
        obj.validMicrosaccadeAmplitudeIndices = setdiff(obj.validMicrosaccadeAmplitudeIndices, idx);
    
        % Compute the displacement vector
        microSaccadeJumpArcMin = amplitudeArcMin*[cos(theta);sin(theta)];
    else
        error('Unknown microsaccade type: ''%s''\n', obj.microSaccadeType);
    end
    
    % Saccade jump in native units
    saccadeJump = microSaccadeJumpArcMin/obj.scalarToArcMin;
    obj.microSaccadeResidualPath = ones(2,durationMilliSeconds);
    obj.microSaccadeResidualPath(1,:) = saccadeJump(1,1)/durationMilliSeconds;
    obj.microSaccadeResidualPath(2,:) = saccadeJump(2,1)/durationMilliSeconds;
    
    % Insert the saccade jump
    endingPosition = obj.emPosTimeSeries(:,tStep+1) + saccadeJump;
    obj.emPosTimeSeries(:,tStep+1) = obj.emPosTimeSeries(:,tStep+1) + obj.microSaccadeResidualPath(:,1);
    % Update the microSaccadeResidualPath
    obj.updateMicroSaccadeResidualPath();
  
    % Keep a record of the saccade onset time indices, saccade jumps and intervals 
    obj.microSaccadeOnsetStepIndices(numel(obj.microSaccadeOnsetStepIndices)+1) = tStep+1;
    saccadeIndex = numel(obj.microSaccadeOnsetStepIndices);
    obj.microSaccadeAmplitudesArcMin(saccadeIndex) = norm(microSaccadeJumpArcMin);
    obj.microSaccadeSpeedsDegsPerSecond(saccadeIndex) = (norm(microSaccadeJumpArcMin)/60)/(durationMilliSeconds/1000);
    obj.microSaccadeIntervalsMilliseconds(saccadeIndex) = intervalSinceLastMicroSaccadeSeconds*1000;
    if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
        obj.microSaccadeTargetLikelihoodMapsList(saccadeIndex,:,:) = targetLikelihoodSpatialMap;
    end
    obj.microSaccadePositions(saccadeIndex,:) = obj.scalarToArcMin * (endingPosition - obj.emPosTimeSeries(:,obj.stabilizationStepsNum))';

    % Update time of last saccade
    obj.lastMicroSaccadeTimeStep = obj.microSaccadeOnsetStepIndices(end);
end

function updateMicroSaccadeResidualPath(obj)
    % Update the microSaccadeResidualPath
    obj.microSaccadeResidualPath = circshift(obj.microSaccadeResidualPath,-1,2);
    obj.microSaccadeResidualPath(:,end) = 0;
end

function targetPosArcMin = posOfMaxLikelihood(obj,targetLikelihoodSpatialMap)
    [~,idx] = max(targetLikelihoodSpatialMap(:));
    if (numel(idx) > 1)
        % More than one locations have the same max. Pick one at random
        index2 = randperm(numel(idx));
        idx = idx(index2(1));
    end
    [row,col] = ind2sub(size(targetLikelihoodSpatialMap), idx);
    targetPosArcMin = [obj.heatMapXYsupport(col) obj.heatMapXYsupport(row)];
end

function currentPositionHeatLevel = updateHeatMap(obj, tStep) 
    % Re-center and scale to ArcMins
    postStabilizationStep = 1+tStep-obj.stabilizationStepsNum;
    xPos = obj.scalarToArcMin*(obj.emPosTimeSeries(1,tStep) - obj.emPosTimeSeries(1,obj.stabilizationStepsNum));
    yPos = obj.scalarToArcMin*(obj.emPosTimeSeries(2,tStep) - obj.emPosTimeSeries(2,obj.stabilizationStepsNum));

    % Compute visited grid coords
    gridNodesNum = length(obj.heatMapXYsupport);
    node0 = (gridNodesNum-1)/2+1;
    visitedGridCol = node0 + sign(xPos) * round(abs(xPos)/obj.heatMapSpatialSampleArcMin);
    visitedGridRow = node0 + sign(yPos) * round(abs(yPos)/obj.heatMapSpatialSampleArcMin);

    % At the beginning of the heatMapUpdateInterval, zero the spatial accumulation map
    if (mod(postStabilizationStep-1,obj.heatMapUpdateIntervalStepsNum) == 0)
        obj.emAccumMap = zeros(gridNodesNum, gridNodesNum);
    end

    % Accumulate the spatial map according to point visited
    if ((visitedGridCol>0) && (visitedGridCol <= gridNodesNum) && (visitedGridRow>0) && (visitedGridRow <= gridNodesNum))
        obj.emAccumMap(visitedGridRow, visitedGridCol) = obj.emAccumMap(visitedGridRow, visitedGridCol)+1;
    end
    
    % At the end of the heatMapUpdateInterval, perform spatial
    % convolution of the accum map with the spatial kernel
    if (mod(postStabilizationStep-1,obj.heatMapUpdateIntervalStepsNum) == obj.heatMapUpdateIntervalStepsNum-1)  
        % Perform spatial convolution with heatMapSpatialKernel
        subSampledPostStabilizationStep = floor((postStabilizationStep-1)/obj.heatMapUpdateIntervalStepsNum)+1;        
        obj.heatMapTimeSeriesIntermediate(subSampledPostStabilizationStep,:,:) = conv2(obj.emAccumMap, obj.heatMapSpatialKernel, 'same');
        
        % Perform temporal convolution with heatMapTemporalKernel
        tBinsForTemporalConvolution = 1:subSampledPostStabilizationStep;
        firstBinAffected = max([1 tBinsForTemporalConvolution(end)-size(obj.heatMapTemporalKernel,1)+1]);
        binsAffected = firstBinAffected : tBinsForTemporalConvolution(end);
        if (numel(binsAffected) == size(obj.heatMapTemporalKernel,1))
            % Use the entire temporal filter
            instantaneousHeatMap = dot(obj.heatMapTimeSeriesIntermediate(binsAffected,:,:), obj.heatMapTemporalKernel, 1);
        else
            % Use the valid portion of the temporal filter
            kernelPortion = obj.heatMapTemporalKernel(1:numel(binsAffected),:,:);
            instantaneousHeatMap = dot(obj.heatMapTimeSeriesIntermediate(binsAffected,:,:), kernelPortion, 1);
        end
        % Inject copies of the instantaneousHeatMap map at appropriate delays into the future
        idx = postStabilizationStep + (0:obj.heatMapUpdateIntervalStepsNum-1);
        obj.heatMapTimeSeries(idx,:,:) = repmat(instantaneousHeatMap, [numel(idx) 1 1]);
    end
    
    % Return heat level at current em position
    currentPositionHeatLevel = 0;
    if ((visitedGridCol>0) && (visitedGridCol <= gridNodesNum) && ...
        (visitedGridRow>0) && (visitedGridRow <= gridNodesNum))
            currentPositionHeatLevel = obj.heatMapTimeSeries(postStabilizationStep,visitedGridRow, visitedGridCol);
    end
    % Keep a history
    obj.currentPositionHeatLevelTimeSeries(1,postStabilizationStep) = currentPositionHeatLevel;
end

function trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds)
    % Trim: only keep samples after the stabilization time
    keptSteps = obj.stabilizationStepsNum+1:obj.tStepsNum-1;
    obj.timeAxis = obj.timeAxis(keptSteps);
    obj.emPosTimeSeries = obj.emPosTimeSeries(:,keptSteps);
    obj.heatMapTimeSeries = obj.heatMapTimeSeries(1:numel(keptSteps),:,:);
    obj.velocityTimeSeries = obj.velocityTimeSeries(keptSteps);
    obj.microSaccadeOnsetStepIndices = obj.microSaccadeOnsetStepIndices-obj.stabilizationStepsNum;
    obj.microSaccadeOnsetStepIndices = obj.microSaccadeOnsetStepIndices(obj.microSaccadeOnsetStepIndices>=0);
    
    % Re-center
    obj.timeAxis = obj.timeAxis - obj.timeAxis(1);
    obj.emPosTimeSeries(1,:) = obj.emPosTimeSeries(1,:) - obj.emPosTimeSeries(1,1);
    obj.emPosTimeSeries(2,:) = obj.emPosTimeSeries(2,:) - obj.emPosTimeSeries(2,1);

    % Resample in time according to passed sampleDurationSeconds
    if (abs(obj.timeStepDurationSeconds - sampleDurationSeconds) > 100*eps(sampleDurationSeconds))
        % Resampled time axis
        oldTimeAxis = obj.timeAxis;
        obj.timeAxis = 0:sampleDurationSeconds:obj.timeAxis(end);
        obj.emPosTimeSeries = obj.smartInterpolation(oldTimeAxis, obj.emPosTimeSeries, obj.timeAxis);
        obj.velocityTimeSeries = obj.smartInterpolation(oldTimeAxis, obj.velocityTimeSeries, obj.timeAxis);
        obj.heatMapTimeSeries = obj.smartInterpolation(oldTimeAxis, obj.heatMapTimeSeries, obj.timeAxis);
        %obj.emPosTimeSeries = (interp1(oldTimeAxis, obj.emPosTimeSeries', obj.timeAxis))';
        %obj.velocityTimeSeries = interp1(oldTimeAxis, obj.velocityTimeSeries, obj.timeAxis);
        %[T,N,M] = size(obj.heatMapTimeSeries);
        %tmp = reshape(obj.heatMapTimeSeries,[T N*M]);
        %obj.heatMapTimeSeries = reshape((interp1(oldTimeAxis, tmp, obj.timeAxis)), [numel(obj.timeAxis) N M]);
    end
end

function outMatrix = smartInterpolation(obj,inputTimeAxis, inputMatrix, outputTimeAxis)
    if (ndims(inputMatrix) > 3)
        error('input matrix must be 1D or 2D')
    elseif (ndims(inputMatrix)==3)
        if (length(inputTimeAxis) == size(inputMatrix,1))
            [T,N,M] = size(obj.heatMapTimeSeries);
            inputMatrix = reshape(inputMatrix,[T N*M]);
            outMatrix = reshape((interp1(inputTimeAxis, inputMatrix, outputTimeAxis)), [numel(outputTimeAxis) N M]);
        elseif (length(inputTimeAxis) == size(inputMatrix,3))
            [N,M,T] = size(obj.heatMapTimeSeries);
            inputMatrix = reshape(inputMatrix,[N*M T]);
            outMatrix = reshape((interp1(inputTimeAxis, inputMatrix', outputTimeAxis))', [numel(outputTimeAxis) N M]);
        else
            error('Time dimension should be either first or last')
        end
    else
        transposeOutMatrix = false;
        if (length(inputTimeAxis) == size(inputMatrix,2))
            inputMatrix = inputMatrix';
            transposeOutMatrix = true;
        end
        outMatrix = interp1(inputTimeAxis, inputMatrix, outputTimeAxis);
        if (transposeOutMatrix)
           outMatrix  = outMatrix'; 
        end
    end
end

% Method to compute the spatial support for the heatMap
function heatMapXYsupport = spatialSupportForHeatMap(obj)
    gridNodesNum = ceil(obj.heatMapWidthArcMin / obj.heatMapSpatialSampleArcMin);
    if (mod(gridNodesNum,2) == 0)
        gridNodesNum = gridNodesNum + 1;
    end
    node0 = (gridNodesNum-1)/2+1;
    heatMapXYsupport = (-(node0-1):(node0-1))*obj.heatMapSpatialSampleArcMin;
end

% Method to compute the spatial kernel for the heat map
function [kernel, kernelSupport] = spatialKernelForHeatMap(obj)
    if (obj.heatMapKernelSpaceConstantArcMin < 0.01)
        sigma = 0.01;
    else
        sigma = obj.heatMapKernelSpaceConstantArcMin;
    end
    halfPixelsNum = round(2.5*sigma/obj.heatMapSpatialSampleArcMin);
    ii = -halfPixelsNum:halfPixelsNum;
    pixels = numel(ii);
    kernel = fspecial('gaussian', pixels, sigma);
    kernelSupport = ii*obj.heatMapSpatialSampleArcMin;
end


% Method to compute the temporal kernel for the heat map
function [kernel, kernelSupport] = temporalKernelForHeatMap(obj, gridNodesNum)
    kernelLength = round(obj.heatMapKernelTimeConstantSeconds/(obj.heatMapUpdateIntervalStepsNum*obj.timeStepDurationSeconds));
    kernelSupport = (0:(kernelLength-1))*obj.heatMapUpdateIntervalStepsNum*obj.timeStepDurationSeconds;
    if (obj.heatMapKernelTimeConstantSeconds < obj.timeStepDurationSeconds)
        kernelTimeConstantSeconds = obj.timeStepDurationSeconds;
    else
        kernelTimeConstantSeconds = obj.heatMapKernelTimeConstantSeconds;
    end
    kernel = exp(-(kernelSupport/kernelTimeConstantSeconds));
    
    % subsample according to obj.heatMapTemporalSampleSeconds
    indices = 1:round(kernelSupport(end)/obj.heatMapTemporalSampleSeconds);
    kernel = kernel(indices);
    kernelSupport = kernelSupport(indices);
    % zero kernel below 1%
    kernel = kernel(kernel>=0.01);
    kernel = (kernel-min(kernel))/(max(kernel)-min(kernel));
    % flip kernel because instead of convolution we are doing a dot product
    kernel = fliplr(kernel(indices));
    kernel = kernel';
    kernelSupport = kernelSupport';
    kernel = repmat(kernel, [1  gridNodesNum gridNodesNum]);
end

% Method to compute obj.microSaccadeTargetLikelihoodSpatialMap
function microSaccadeTargetLikelihoodSpatialMap = computeMicrosaccadeTargetLikelihoodSpatialMap(obj, heatMapXYsupportArcMin)
    % Gaussian with sigma obj.microaccadeTargetSigma
    [x,y] = meshgrid(heatMapXYsupportArcMin,heatMapXYsupportArcMin);
    microSaccadeTargetLikelihoodSpatialMap = exp(-0.5*(x/obj.fixationMapSpaceConstantArcMin).^2) .* exp(-0.5*(y/obj.fixationMapSpaceConstantArcMin).^2);
    microSaccadeTargetLikelihoodSpatialMap = microSaccadeTargetLikelihoodSpatialMap / max(microSaccadeTargetLikelihoodSpatialMap(:));
end

end % private methods

end % class