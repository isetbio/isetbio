function compute(obj, emDurationSeconds, sampleDurationSeconds, ...
    nTrials, computeVelocity, varargin)
% Compute nTrials of fixational eye movements.
%
% Syntax:
%   compute(obj, emDurationSeconds, sampleDurationSeconds, ...
%       nTrials, computeVelocity, [varargin])
%
% Description:
%    Compute the fixational eye movements for a fixationalEM object.
%
% Inputs:
%    obj                   - Object. The fixationalEM object.
%    emDurationSeconds     - Numeric. Eye movement duration, in seconds.
%    sampleDurationSeconds - Numeric. The duration of the sample period,
%                            also in seconds.
%    nTrials               - Numeric. The number of trials.
%    centerPaths           - Logical. Whether to center the eye movement paths
%    computeVelocity       - Boolean. A boolean indicating whether or not
%                            to also compute the velocity.
%    varargin              - (Optional) Additional parameter(s) that may be
%                            required to execute the function.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    'useParFor'           - Boolean. Whether to compute using a parfor 
%                            loop (to take advantage of multiple processors)
%                            of to use a single processor. Default false,
%                            i.e., use a single processor.
%
% For usage see:
%  t_fixationalEyeMovementsTypes
%
% History:
%    01/03/18  NPC      ISETBIO Team, 2018
%    05/15/18  jnm      Formatting
%    05/24/18  BW, NPC  Comments
%    6/04/19   NPC      Added centerPaths key/value pair

% Parse inputs
p = inputParser;
p.addRequired('emDurationSeconds', @isnumeric);
p.addRequired('sampleDurationSeconds', @isnumeric);
p.addRequired('nTrials', @isnumeric);
p.addRequired('computeVelocity', @islogical);
p.addParameter('useParfor', false, @islogical);
p.addParameter('centerPaths', false, @islogical);
p.parse(emDurationSeconds, sampleDurationSeconds, ...
    nTrials, computeVelocity, varargin{:});

% Reset output arrays
obj.initOutputs();

% Set random seed
if (isempty(obj.randomSeed))
    %rng('shuffle');                         % 2018b and before
    rng('shuffle', 'CombRecursive');         % 2019a and after
else
    %rng(obj.randomSeed);                    % 2018b and before
    rng(obj.randomSeed, 'CombRecursive');    % 2019a and after
end

if (obj.beVerbose)
    fprintf(['Computing emModel for %2.2f seconds; sampleT = ' ...
        '%2.4f sec, samples: %d\n'], emDurationSeconds, ...
        sampleDurationSeconds, obj.tStepsNum);
end

% Compute first trial
iTrial = 1;
computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds);
nTimeBins = length(obj.timeAxis);

allTrialsEmPosArcMin = zeros(nTrials, nTimeBins, 2);
allTrialsEmPosArcMin(iTrial, :, :) = ...
    reshape(obj.emPosTimeSeriesArcMin', ...
    [1 nTimeBins 2]);


if (computeVelocity)
    allTrialsVelocityArcMin = zeros(nTrials, ...
        length(obj.velocityArcMinPerSecTimeSeries), 1);
    allTrialsVelocityArcMin(iTrial, :) = ...
        obj.velocityArcMinPerSecTimeSeries;
end

% Compute remaining trials
if (p.Results.useParfor)
    parfor iTrial = 2:nTrials
        computeSingleTrial(...
            obj, emDurationSeconds, sampleDurationSeconds);
        allTrialsEmPosArcMin(iTrial, :, :) = ...
            reshape(obj.emPosTimeSeriesArcMin', ...
            [1 length(obj.emPosTimeSeriesArcMin) 2]);
        if (computeVelocity)
            allTrialsVelocityArcMin(iTrial, :) = ...
                obj.velocityArcMinPerSecTimeSeries;
        end
    end
else
    for iTrial = 2:nTrials
        computeSingleTrial(...
            obj, emDurationSeconds, sampleDurationSeconds);
        allTrialsEmPosArcMin(iTrial, :, :) = ...
            reshape(obj.emPosTimeSeriesArcMin', ...
            [1 length(obj.emPosTimeSeriesArcMin) 2]);
        if (computeVelocity)
            allTrialsVelocityArcMin(iTrial, :) = ...
                obj.velocityArcMinPerSecTimeSeries;
        end
    end
end

% Center all paths on the origin
if (p.Results.centerPaths)
    centers = mean(allTrialsEmPosArcMin,2);
    allTrialsEmPosArcMin = bsxfun(@minus, allTrialsEmPosArcMin, centers);
end

obj.emPosArcMin = allTrialsEmPosArcMin;
if (computeVelocity), obj.velocityArcMin = allTrialsVelocityArcMin; end

end