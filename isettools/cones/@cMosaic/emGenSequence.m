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
    p.addParameter('nTrials', 1, @isscalar);
    p.addParameter('useParfor', false, @islogical);
    p.addParameter('randomSeed', [], @isscalar);
    p.parse(obj, durationSeconds, varargin{:});
    
    % Save current random seed
    oldRandomSeed = rng;
    
    % Compute number of eye movements
    nEyeMovements = round(p.Results.durationSeconds / obj.integrationTime);
    
    % Instantiate fixationaEM object and attach it to the cMosaic object
    obj.fixEMobj = fixationalEM();
    obj.fixEMobj.microSaccadeType = p.Results.microsaccadeType;
    if (~isempty(p.Results.microsaccadeMeanIntervalSeconds))
        obj.fixEMobj.microSaccadeMeanIntervalSeconds = p.Results.microsaccadeMeanIntervalSeconds;
    end
    obj.fixEMobj.computeForCmosaic(obj, nEyeMovements, ...
        'nTrials', p.Results.nTrials, ...
        'rSeed', p.Results.randomSeed, ...
        'useParfor', p.Results.useParfor);
    
    % Restore previous random seed
    rng(oldRandomSeed);
    
end
