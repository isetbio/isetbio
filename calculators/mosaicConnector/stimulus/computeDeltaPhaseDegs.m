function [deltaPhaseDegs, actualTemporalFrequencyHz] = computeDeltaPhaseDegs(temporalFrequencyHz, mosaicIntegrationTimeSeconds)
    mosaicIntegrationTimeMilliseconds = mosaicIntegrationTimeSeconds * 1000;
    periodMilliseconds = 1000/temporalFrequencyHz;
    [scenesNum, actualTemporalFrequencyHz] = scenesNumForPeriod(periodMilliseconds, mosaicIntegrationTimeMilliseconds);
    k = 1;
    while (scenesNum > 20) && (k < 5)
        k = k + 1;
        % see if we can get by using double the integrationTime
        [scenesNum, actualTemporalFrequencyHz] = scenesNumForPeriod(periodMilliseconds, k*mosaicIntegrationTimeMilliseconds);
    end
    deltaPhaseDegs = 360/scenesNum;
    frameDurationMilliseconds = periodMilliseconds/scenesNum;
end

function [scenesNum, actualTemporalFrequencyHz] = scenesNumForPeriod(periodMilliseconds, mosaicIntegrationTimeMilliseconds)
    periodMilliseconds = round(periodMilliseconds/mosaicIntegrationTimeMilliseconds)*mosaicIntegrationTimeMilliseconds;
    actualTemporalFrequencyHz = 1000.0/periodMilliseconds;
    scenesNum = periodMilliseconds/mosaicIntegrationTimeMilliseconds;
end
