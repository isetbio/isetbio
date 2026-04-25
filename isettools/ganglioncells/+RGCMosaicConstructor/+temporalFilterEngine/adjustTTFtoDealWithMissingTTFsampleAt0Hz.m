%
% RGCMosaicConstructor.temporalFilterEngine.adjustTTFtoDealWithMissingTTFsampleAt0Hz()
%

function [theTTF,  temporalFrequenciesExamined] = adjustTTFtoDealWithMissingTTFsampleAt0Hz(...
    theTTF,  temporalFrequenciesExamined, minimumDelaySecondsForEstimationOfBaseline)

    if (temporalFrequenciesExamined(1) == 0)
        fprintf(2, 'There is a 0 Hz point in the temporal frequency support. No adjustment of TTF is needed. \n');
        return;
    end

    % We do not measure a response at 0 Hz. But we need to put something at
    % 0 Hz. We start by adding a point whose amplitude is the same as the
    % amplitude at the first non-zero frequency.

    % Ensure that the distance from the 1st point to the 0Hz point is the
    % same as the frequency sampling interval.
    deltaFrequency = temporalFrequenciesExamined(2)-temporalFrequenciesExamined(1);
    if (temporalFrequenciesExamined(1)-deltaFrequency ~= 0)
        temporalFrequenciesExamined(1)-deltaFrequency
        error('Cannot add 0 TF')
    end

    % Assume TTF(0) = norm(TTF(1))
    temporalFrequenciesExamined = cat(2, 0, temporalFrequenciesExamined);
    theTTF = cat(2, norm(theTTF(1)), theTTF);


    % This introduces a DC term. Estimate it from the impulse response
    thePhotocurrentBasedImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theTTF , temporalFrequenciesExamined, ...
                'window', 'hamming', ...
                'causal', false);

    % We estimate the baseline from time bins > 0.8 minimumDelaySecondsForEstimationOfBaselin
    idx = find(abs(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds)>=minimumDelaySecondsForEstimationOfBaseline);
    baselineOffset = mean(thePhotocurrentBasedImpulseResponseData.amplitude(idx));

    % Correction for the thePhotocurrentBasedMRGCcellTTF(1) 
    dcCorrection = baselineOffset * numel(theTTF)*2;
    theTTF(1) = theTTF(1) - dcCorrection;
end

