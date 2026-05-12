%
% RGCMosaicConstructor.temporalFilterEngine.adjustTTFtoDealWithMissingTTFsampleAt0Hz()
%

function [theTTF,  temporalFrequenciesExamined] = adjustTTFtoDealWithMissingTTFsampleAt0Hz(...
    theTTF,  temporalFrequenciesExamined)

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
    theDCterm = norm(theTTF(1));

    % Or assume TTF(0) = 0;
    theDCterm = 0;

    temporalFrequenciesExamined = cat(2, 0, temporalFrequenciesExamined);
    theTTF = cat(2, theDCterm, theTTF);
end

