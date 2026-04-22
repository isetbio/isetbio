%
% RGCMosaicConstructor.temporalFilterEngine.impulseResponseFunctionFromTTF
%

function theImpulseResponseFunctionStruct = impulseResponseFunctionFromTTF(...
    theTTF, temporalFrequencySupportHz, performFFTshift, zeroPaddingLength)


    if (temporalFrequencySupportHz(1) > 0)
        % Add the zero TF point
        deltaFrequency = temporalFrequencySupportHz(2)-temporalFrequencySupportHz(1);
        if (temporalFrequencySupportHz(1)-deltaFrequency ~= 0)
            temporalFrequencySupportHz(1)-deltaFrequency
            error('Cannot add 0 TF')
        end
        % Assume TTF(0) = TTF(1st TF bin)
        temporalFrequencySupportHz = cat(2, 0, temporalFrequencySupportHz);
        theTTF = cat(2, theTTF(1), theTTF);
    end

    % Zero padding
    if (~isempty(zeroPaddingLength))
        extraSamplesNum = zeroPaddingLength-numel(temporalFrequencySupportHz);
        nSamples = numel(theTTF) + extraSamplesNum;
        theTTF = cat(2, theTTF, zeros(1,nSamples-numel(theTTF)));
        temporalFrequencySupportHz = temporalFrequencySupportHz(1) + (0:(nSamples-1))*(temporalFrequencySupportHz(2)-temporalFrequencySupportHz(1));
    end

    
    % Convert single-sided spectrum to double sided
    theDoubleSidedTTF = [theTTF(1) theTTF(2:end)/2 fliplr(conj(theTTF(2:end)))/2];

    % Inverse FFT
    theImpulseResponse = ifft(theDoubleSidedTTF, 'symmetric');
    if (performFFTshift)
        theImpulseResponse = fftshift(theImpulseResponse);
    end

    % Nyquist frequency
    fMax = max(temporalFrequencySupportHz);
    dtSeconds = 1/(2*fMax);
    theTemporalSupportSeconds = (0:(numel(theImpulseResponse)-1)) * dtSeconds;

    theImpulseResponseFunctionStruct.amplitude = theImpulseResponse;
    theImpulseResponseFunctionStruct.temporalSupportSeconds = theTemporalSupportSeconds;
end