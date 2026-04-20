%
% RGCMosaicConstructor.temporalFilterEngine.impulseResponseFunctionFromTTF
%

function theImpulseResponseFunctionStruct = impulseResponseFunctionFromTTF(theTTF, temporalFrequencySupportHz)

    % Zero padding
    extraSamplesNum = 0;
    nSamples = numel(theTTF) + extraSamplesNum;
    theTTF = cat(2, theTTF, zeros(1,nSamples-numel(theTTF)));
    temporalFrequencySupportHz = temporalFrequencySupportHz(1) + (0:(nSamples-1))*(temporalFrequencySupportHz(2)-temporalFrequencySupportHz(1));

    % Convert single-sided spectrum to double sided
    theDoubleSidedTTF = [theTTF(1) theTTF(2:end)/2 fliplr(conj(theTTF(2:end)))/2];

    % Inverse FFT
    theImpulseResponse = ifft(theDoubleSidedTTF, 'symmetric');

    % Nyquist frequency
    fMax = max(temporalFrequencySupportHz);
    dtSeconds = 1/(2*fMax);
    theTemporalSupportSeconds = (0:(numel(theImpulseResponse)-1)) * dtSeconds;

    theImpulseResponseFunctionStruct.amplitude = theImpulseResponse;
    theImpulseResponseFunctionStruct.temporalSupportSeconds = theTemporalSupportSeconds;
end