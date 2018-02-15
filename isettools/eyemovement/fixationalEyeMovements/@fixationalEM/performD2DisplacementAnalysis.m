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