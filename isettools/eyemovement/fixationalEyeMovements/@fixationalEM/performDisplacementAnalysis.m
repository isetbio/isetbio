function [meanD, maxD, meanDscrambled, DrandomWalk, ...
    timeLagsMilliseconds] = performDisplacementAnalysis(emPos, ...
    timeAxisSeconds, varargin)
% Run the displacement analysis
%
% Syntax:
%   [meanD, maxD, meanDscrambled, DrandomWalk, timeLagMilliseconds] = ...
%       performDisplacementAnalysis(emPos, timeAxisSeconds, [varargin])
%
% Description:
%    Perform the displacement analysis using the provided information.
%
%    This code contains examples. To access, type 'edit
%    performDisplacementAnalysis.m' into the Command Window.
%
% Inputs:
%    emPos               - Vector. Eye movement positions.
%    timeAxisSeconds     - Vector. The time axis vector.
%    varargin            - (Optional) Additional argument(s) as needed, to
%                          perform the analysis. Usually contained in
%                          key/value pairs (see section below).
%
% Outputs:
%    meanD               - Vector. The mean displacements.
%    maxD                - Vector. The maximum displacements.
%    meanDscrambled      - Vector. The scrambled mean displacements.
%    DrandomWalk         - Vector. Random walk of displacements.
%    timeLagMilliseconds - Vector. The time lag for displacements in ms.
%
% Optional key/value pairs:
%    'mode'              - String. A string describing the mode. Options
%                          are 'D1', and 'D2'. Default 'D1'.
%

% Examples:
%{
    emDurationSeconds = 2;
    sampleTimeSeconds = .001;
    nTrials = 512;

    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = 'stats based';

    % Compute the emPaths
    computeVelocity = true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
        nTrials, computeVelocity, 'useParfor', true);

    emPathArcMin = squeeze(fixEMobj.emPosArcMin(1,:,:));
    xPosArcMin = squeeze(emPathArcMin(:,1));

    [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        fixEMobj.performDisplacementAnalysis(xPosArcMin', ...
        fixEMobj.timeAxis, 'mode', 'D1');
%}
p = inputParser;
addParameter(p, 'mode', 'D1', ...
    @(x) (ischar(x) && ismember(x, {'D1', 'D2'})));
parse(p, varargin{:});

timeLagsMilliseconds = [1 2:2:100 105:5:1000];
dtMilliseconds = (timeAxisSeconds(2) - timeAxisSeconds(1)) * 1000;

% Compute emPos from randomized intervals
sdiff = diff(emPos);
scrambledIndices = randperm(length(sdiff));
sdiff = sdiff(scrambledIndices);
emPosScrambled = [emPos(1) emPos(1) + cumsum(sdiff)];

% Initialize displacement vectors
meanD = zeros(1, numel(timeLagsMilliseconds));
maxD = zeros(1, numel(timeLagsMilliseconds));
meanDscrambled = zeros(1, numel(timeLagsMilliseconds));

tStepsNum = length(timeAxisSeconds);
for delayIndex = 1:numel(timeLagsMilliseconds)
    m = round(timeLagsMilliseconds(delayIndex) / dtMilliseconds);
    sumDeltaRsquared = zeros(1, tStepsNum - m);
    sumDeltaRsquaredScrambled = zeros(1, tStepsNum - m);
    for i = 1:tStepsNum - m
        sumDeltaRsquared(i) = abs(emPos(i + m) - emPos(i));
        sumDeltaRsquaredScrambled(i) = abs(emPosScrambled(i + m) - ...
            emPosScrambled(i));
    end
    if (strcmp(p.Results.mode, 'D2'))
        sumDeltaRsquared = sumDeltaRsquared .^ 2;
        sumDeltaRsquaredScrambled = sumDeltaRsquaredScrambled .^ 2;
    end
    if (numel(sumDeltaRsquared) > 0)
        maxD(1, delayIndex) = max(sumDeltaRsquared);
        meanD(1, delayIndex) = mean(sumDeltaRsquared);
        meanDscrambled(1, delayIndex) = mean(sumDeltaRsquaredScrambled);
    end
end
k = meanDscrambled(1) / timeLagsMilliseconds(1);
H = 0.5;
DrandomWalk = k * timeLagsMilliseconds .^ H;

end