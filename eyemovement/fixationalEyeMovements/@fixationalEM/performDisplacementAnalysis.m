function [meanD, maxD, meanDscrambled, DrandomWalk, ...
    timeLagsMilliseconds] = performDisplacementAnalysis(emPos, ...
    timeAxisSeconds, varargin)
% Run the displacement analysis for passed emPos component (x or y)
%
% Syntax:
%   [meanD, maxD, meanDscrambled, DrandomWalk, timeLagMilliseconds] = ...
%       performDisplacementAnalysis(emPos, timeAxisSeconds, [varargin])
%
% Description:
%    Perform the D1 or the D2 displacement analysis using the provided 
%    information.
%
%    This code contains examples. To access, type 'edit
%    performDisplacementAnalysis.m' into the Command Window.
%
% Inputs:
%    emPos               - Vector. Eye movement positions.
%    timeAxisSeconds     - Vector. The time axis of the em positions
%    varargin            - (Optional) Additional argument(s) as needed, to
%                          perform the analysis. Usually contained in
%                          key/value pairs (see section below).
%
% Outputs:
%    meanD               - Vector. Mean displacement function for the
%                          passed emPath.
%    maxD                - Vector. Max displacement function for the 
%                          passed emPath.
%    meanDscrambled      - Vector. Mean displacement function for the 
%                          path constructed by scrambling the intervals.
%    DrandomWalk         - Vector. Mean isplacement function for a 
%                          random walk emPath.
%    timeLagMilliseconds - Vector. The support of the displacement function
%
% Optional key/value pairs:
%    'mode'              - String. A string describing the mode. Options
%                          are 'D1', and 'D2'. Default 'D1'.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments


% Examples:
%{
    emDurationSeconds = 2;
    sampleTimeSeconds = .001;
    nTrials = 512;

    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = 'none';

    % Compute the emPaths
    computeVelocity = true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
        nTrials, computeVelocity, 'useParfor', true);

    for trialNo = 1:nTrials
        xPosArcMin = squeeze(fixEMobj.emPosArcMin(trialNo,:,1));
        [meanD2(trialNo,:), ~, meanD2scrambled(trialNo,:), ...
         D2randomWalk, timeLagsMilliseconds] = ...
            fixEMobj.performDisplacementAnalysis(xPosArcMin, ...
        fixEMobj.timeAxis, 'mode', 'D2');
    end

    meanD2 = mean(meanD2,1);
    meanD2scrambled = mean(meanD2scrambled,1);

    figure();
    plot(timeLagsMilliseconds, meanD2, 'k-'); hold on;
    plot(timeLagsMilliseconds, meanD2scrambled, 'k--'); hold off
    set(gca, 'XScale', 'log', 'Yscale', 'log'); axis 'square';
    xlabel('interval (ms)'); ylabel('mean D2 displacement (arc min^2)');
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