function [meanD, maxD, meanDscrambled, DrandomWalk, ...
    timeLagsMilliseconds] = performD2DisplacementAnalysis(...
    emPos, timeAxisSeconds)
% Run the delta squared displacement analysis
%
% Syntax:
%   [meanD, maxD, meanDscrambled, DrandomWalk, timeLagMilliseconds] = ...
%       performD2DisplacementAnalysis(emPos, timeAxisSeconds)
%
% Description:
%    Perform the delta squared displacement analysis.
%
%    This code contains examples. To access, type 'edit
%    performDisplacementAnalysis.m' into the Command Window.
%
% Inputs:
%    emPos               - Vector. Eye movement positions.
%    timeAxisSeconds     - Vector. The time axis vector in seconds.
%
% Outputs:
%    meanD               - Vector. The mean d^2 displacements.
%    maxD                - Vector. The maximum d^2 displacements.
%    meanDscrambled      - Vector. The scrambled mean d^2 displacements.
%    DrandomWalk         - Vector. Random walk of d^2 displacements.
%    timeLagMilliseconds - Vector. The time lag for d^2 displacements in ms
%
% Optional key/value pairs:
%    None.
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

    emPathArcMin = squeeze(fixEMobj.emPosArcMin(1, :, :));
    xPosArcMin = squeeze(emPathArcMin(:, 1));

    [meanD, maxD, meanDscrambled, DrandomWalk, timeLagsMilliseconds] = ...
        fixEMobj.performD2DisplacementAnalysis(xPosArcMin', ...
        fixEMobj.timeAxis);
%}

timeLagsMilliseconds = [1 2:2:100 105:5:1000];
dtMilliseconds = (timeAxisSeconds(2) - timeAxisSeconds(1)) * 1000;

% Compute emPos from randomized intervals
sdiff = diff(emPos);
scrambledIndices = randperm(length(sdiff));
sdiff = sdiff(scrambledIndices);
emPosScrambled = [emPos(1) emPos(1) + cumsum(sdiff)];

% Initialize displacement vectors
meanD = nan(1, numel(timeLagsMilliseconds));
maxD = nan(1, numel(timeLagsMilliseconds));
meanDscrambled = nan(1, numel(timeLagsMilliseconds));

tStepsNum = length(timeAxisSeconds);
for delayIndex = 1:numel(timeLagsMilliseconds)
    m = round(timeLagsMilliseconds(delayIndex) / dtMilliseconds);
    sumDeltaRsquared = zeros(1, tStepsNum - m);
    sumDeltaRsquaredScrambled = zeros(1, tStepsNum - m);
    for i = 1:tStepsNum-m
        sumDeltaRsquared(i) = (norm(emPos(i + m) - emPos(i))) ^ 2;
        sumDeltaRsquaredScrambled(i) = ...
            (norm(emPosScrambled(i + m) - emPosScrambled(i))) ^ 2;
    end
    if (numel(sumDeltaRsquared) > 0)
        maxD(delayIndex) = max(sumDeltaRsquared);
        meanD(delayIndex) = mean(sumDeltaRsquared);
        meanDscrambled(delayIndex) = mean(sumDeltaRsquaredScrambled);
    end
end
k = meanDscrambled(1) / timeLagsMilliseconds(1);
H = 0.5;
DrandomWalk = k * timeLagsMilliseconds .^ H;

end