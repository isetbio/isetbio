%
% RGCMosaicConstructor.temporalFilterEngine.oneSidedTTFfromTemporalImpulseResponse
%

function theOneSidedTTF = oneSidedTTFfromTemporalImpulseResponse(theIR)

    theTwoSidedTTF = fft(theIR);

    N = numel(theTwoSidedTTF);
    M = N/2 + 1;
    % Generate one sided TTF
    theOneSidedTTF(1:M) = theTwoSidedTTF(1:M);
    theOneSidedTTF(2:M) = 2*theTwoSidedTTF(2:M);

end
