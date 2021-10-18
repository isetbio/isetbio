function emGenSequence(obj, durationSeconds, varargin)
% Generate sequence of fixational eye movements for a cMosaic

    % Help
    if (ischar(durationSeconds)) && strcmp(durationSeconds,'help')
        doc('cMosaic.emGenSequence');
        return;
    end

    % Parse input
    p = inputParser;
    p.addRequired('obj', @(x)(isa(x,'cMosaic')));
    p.addRequired('durationSeconds', @isscalar);
    validMicroSaccadeTypes = {'none','stats based','heatmap/fixation based'};
    p.addParameter('microsaccadeType', 'none', @(x)(ismember(x,validMicroSaccadeTypes)));
    p.addParameter('microsaccadeMeanIntervalSeconds', [], @isscalar);
    p.addParameter('centerPaths', false, @islogical);
    p.addParameter('centerPathsAtSpecificTimeMsec', [], @isnumeric);
    p.addParameter('centerPathsAtSpecificTime', [], @isnumeric);
    p.addParameter('nTrials', 1, @isscalar);
    p.addParameter('computeVelocity', false, @islogical);
    p.addParameter('useParfor', false, @islogical);
    p.addParameter('randomSeed', [], @isscalar);
    p.addParameter('driftModelControlGain', [], @isscalar);
    p.addParameter('driftModelFeedbackGain', [], @isscalar);
    p.addParameter('driftModelFeedbackXDelaySeconds', [], @isscalar);
    p.addParameter('driftModelFeedbackYDelaySeconds', [], @isscalar);
    p.addParameter('driftModelPositionNoiseStd', [], @isscalar);
    
    p.parse(obj, durationSeconds, varargin{:});
    
    % Save current random seed
    oldRandomSeed = rng;
    
    % Compute number of eye movements
    nEyeMovements = round(p.Results.durationSeconds / obj.integrationTime);
    
    if (nEyeMovements < 1)
        error('Zero eye movements. Either increase the eye movement duration (current: %2.0f msec) or decrease the integration time of the cone mosaic (current: %2.0 msec).\n', ...
            p.Results.durationSeconds*1000, obj.integrationTime*1000);
    end
    
    % Instantiate fixationaEM object and attach it to the cMosaic object
    obj.fixEMobj = fixationalEM();
    

    % Change drift model parameters
    if (~isempty(p.Results.driftModelControlGain))
        fprintf('Changed drift CONTROL GAIN from %f to %f\n', ...
            obj.fixEMobj.controlGamma, p.Results.driftModelControlGain);
        obj.fixEMobj.controlGamma = p.Results.driftModelControlGain;
        
    end
    
    if (~isempty(p.Results.driftModelFeedbackXDelaySeconds))
        fprintf('Changed drift FEEDBACK DELAY (X) from %f to %f\n', ...
            obj.fixEMobj.feedbackXposDelaySeconds, p.Results.driftModelFeedbackXDelaySeconds);
        obj.fixEMobj.feedbackXposDelaySeconds = p.Results.driftModelFeedbackXDelaySeconds;
    end
    
    if (~isempty(p.Results.driftModelFeedbackYDelaySeconds))
        fprintf('Changed drift FEEDBACK DELAY (Y) from %f to %f\n', ...
            obj.fixEMobj.feedbackYposDelaySeconds, p.Results.driftModelFeedbackYDelaySeconds);
        obj.fixEMobj.feedbackYposDelaySeconds = p.Results.driftModelFeedbackYDelaySeconds;
    end
    
    if (~isempty(p.Results.driftModelFeedbackGain))
        fprintf('Changed drift FEEDBACK GAIN from %f to %f\n', ...
            obj.fixEMobj.feedbackGain, p.Results.driftModelFeedbackGain);
        obj.fixEMobj.feedbackGain = p.Results.driftModelFeedbackGain;
    end
    
    if (~isempty(p.Results.driftModelPositionNoiseStd))
        fprintf('Changed drift POSITION NOISE STD from %f to %f\n', ...
            obj.fixEMobj.positionNoiseStd, p.Results.driftModelPositionNoiseStd);
        obj.fixEMobj.positionNoiseStd = p.Results.driftModelPositionNoiseStd;
    end
    
    
    % Change micro-saccade parameters
    obj.fixEMobj.microSaccadeType = p.Results.microsaccadeType;
    if (~isempty(p.Results.microsaccadeMeanIntervalSeconds))
        obj.fixEMobj.microSaccadeMeanIntervalSeconds = p.Results.microsaccadeMeanIntervalSeconds;
    end
    
    obj.fixEMobj.computeForCmosaic(obj, nEyeMovements, ...
        'centerPaths', p.Results.centerPaths, ...
        'centerPathsAtSpecificTimeMsec', p.Results.centerPathsAtSpecificTimeMsec, ...
        'nTrials', p.Results.nTrials, ...
        'rSeed', p.Results.randomSeed, ...
        'computeVelocity', p.Results.computeVelocity, ...
        'useParfor', p.Results.useParfor);
    
    % Restore previous random seed
    rng(oldRandomSeed);
    
end
