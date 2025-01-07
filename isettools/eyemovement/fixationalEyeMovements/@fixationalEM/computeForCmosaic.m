function computeForCmosaic(obj, cMosaicOBJ, eyeMovementsPerTrial, varargin)
% Compute fixational eye movement for a @cMosaic cone mosaic
%
% Syntax:
%   computeForCmosaic(obj, cMosaicOBJ, eyeMovementsPerTrial, [varargin])
%   obj.computeForCmosaic(cMosaicOBJ, eyeMovementsPerTrial, [varargin])
%
% Description:
%    Compute the fixational eye movement for a provided cone mosaic.
%
% Inputs:
%    obj                  - Object. The fixationalEM object.
%    cMosaicOBJ           - Object. A cMosaic object.
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
%    'centerPathsAtSpecificTimeMsec - Numeric. When a non-empty time is passed
%                            the emPaths at that specific time are at (0,0)
%                           Default: empty, i.e. dont center.
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
p.addRequired('cMosaicOBJ', @(x)(isa(x, 'cMosaic')));
p.addRequired('eyeMovementsPerTrial', @isscalar);
p.addParameter('nTrials', 1, @isscalar);
p.addParameter('centerPaths', false, @islogical);
p.addParameter('centerPathsAtSpecificTimeMsec', [], @isnumeric);
p.addParameter('computeVelocity', false, @islogical);
p.addParameter('rSeed', [], @(x) (isempty(x) | isscalar(x)));
p.addParameter('useParfor', false, @islogical);
p.parse(cMosaicOBJ, eyeMovementsPerTrial, varargin{:});

% Set optional parameters based on input
obj.randomSeed = p.Results.rSeed;
nTrials = p.Results.nTrials;
centerPaths = p.Results.centerPaths;
centerPathsAtSpecificTimeMsec = p.Results.centerPathsAtSpecificTimeMsec;
computeVelocitySignal = p.Results.computeVelocity;
useParfor = p.Results.useParfor;

% Extract sampleDuration from cMosaic's integration time
sampleDurationSeconds = cMosaicOBJ.integrationTime;

% Compute emDurationSeconds
emDurationSeconds = eyeMovementsPerTrial * sampleDurationSeconds;

% Generate fixational eye movements (thie generates obj.emPosArcMin)
obj.compute(emDurationSeconds, sampleDurationSeconds, nTrials, ...
    computeVelocitySignal, 'useParfor', useParfor, ...
    'centerPaths', centerPaths, 'centerPathsAtSpecificTimeMsec', centerPathsAtSpecificTimeMsec);

% Also return the path in units of microns
obj.emPosMicrons = cMosaicOBJ.distanceDegreesToDistanceMicronsForCmosaic(obj.emPosArcMin / 60);

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