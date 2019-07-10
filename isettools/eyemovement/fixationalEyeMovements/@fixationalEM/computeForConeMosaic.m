function computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, varargin)
% Compute fixational eye movement for a cone mosaic
%
% Syntax:
%   computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, [varargin])
%   obj.computeForConeMosaic(coneMosaic, eyeMovementsPerTrial, [varargin])
%
% Description:
%    Compute the fixational eye movement for a provided cone mosaic.
%
% Inputs:
%    obj                  - Object. The fixationalEM object.
%    coneMosaic           - Object. A coneMosaic object.
%    eyeMovementsPerTrial - Numeric. The number of eye movements per trial.
%    varargin             - (Optional) Additional arguments as desired, see
%                           the optional key/values section.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    'nTrials'            - Numeric. The number of trials. Default 1.
%    'centerPaths'        - Logical. Whether to center the eye movement paths
%    'computeVelocity'    - Boolean. Whether should calculate the velocity.
%                           Default false.
%    'rSeed'              - Numeric. Random seed. Empty causes rng to be inititialized
%                           with 'shuffle', otherwise with the integer seed set here.
%                           Default empty ([]).
%    'useParFor'          - Boolean. <Unknown>. Default false.
%

% History:
%   6/25/18  dhb      Change arg check for rSeed so that empty is allowable.
%                     Update header commment to explain its behavior.
%   6/04/19  npc      Added centerPaths key/value pair

p = inputParser;
p.addRequired('coneMosaic', @(x)(isa(x, 'coneMosaic')));
p.addRequired('eyeMovementsPerTrial', @isscalar);
p.addParameter('nTrials', 1, @isscalar);
p.addParameter('centerPaths', false, @islogical);
p.addParameter('computeVelocity', false, @islogical);
p.addParameter('rSeed', [], @(x) (isempty(x) | isscalar(x)));
p.addParameter('useParfor', false, @islogical);
p.parse(coneMosaic, eyeMovementsPerTrial, varargin{:});

% Set optional parameters based on input
obj.randomSeed = p.Results.rSeed;
nTrials = p.Results.nTrials;
centerPaths = p.Results.centerPaths;
computeVelocitySignal = p.Results.computeVelocity;
useParfor = p.Results.useParfor;

% Extract sampleDuration from coneMosaic's integration time
sampleDurationSeconds = coneMosaic.integrationTime;

% Compute emDurationSeconds
emDurationSeconds = eyeMovementsPerTrial * sampleDurationSeconds;

% Generate fixational eye movements (obj.emPathsArcMin)
compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials, ...
    computeVelocitySignal, 'centerPaths', centerPaths, 'useParfor', useParfor);

% Subsample to cone positions
conePatternSampleMicrons = coneMosaic.patternSampleSize(1) * 1e6;
emPosMicronsNotRounded = ...
    (obj.emPosArcMin / 60) * coneMosaic.micronsPerDegree;
% Return the paths in units of cone mosaic pattern size
obj.emPos = sign(emPosMicronsNotRounded) .* ...
    round(abs(emPosMicronsNotRounded) / conePatternSampleMicrons);
% Also return the path in units of microns
obj.emPosMicrons = obj.emPos * conePatternSampleMicrons;
obj.emPosArcMin = obj.emPosMicrons / coneMosaic.micronsPerDegree * 60;

if (computeVelocitySignal)
    iTrial = 1;
    velocityArcMinTrial1 = ...
        obj.computeVelocity(squeeze(obj.emPosArcMin(iTrial, :, :)));
    allTrialsVelocityArcMin = ...
        zeros(nTrials, length(velocityArcMinTrial1), 1);
    allTrialsVelocityArcMin(iTrial, :) = velocityArcMinTrial1;

    allTrialsEMPosArcMin = obj.emPosArcMin;
    if (useParfor)
        parfor iTrial = 2:nTrials
            allTrialsVelocityArcMin(iTrial, :) = ...
                fixationalEM.computeVelocity(squeeze(...
                allTrialsEMPosArcMin(iTrial, :, :)));
        end
    else
        for iTrial = 2:nTrials
            allTrialsVelocityArcMin(iTrial, :) = obj.computeVelocity(...
                squeeze(allTrialsEMPosArcMin(iTrial, :, :)));
        end
    end
    obj.velocityArcMin = allTrialsVelocityArcMin;
end

end