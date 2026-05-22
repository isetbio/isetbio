%
% RGCMosaicConstructor.temporalFilterEngine.centerAndWindowTemporalImpulseResponse
%
%
function  [theCenteredWaveform, theDelaySeconds] = centerAndWindowTemporalImpulseResponse(temporalSupportSeconds, theWaveform)

error('Are we using this?')

    % Find the time of its peak
    [~, tBinOfMaxResponse] = max(abs(theWaveform(:)));
    
    % Center the impulse response on the time window
    midPoint = round(numel(temporalSupportSeconds)/2);
    theDelaySeconds = temporalSupportSeconds(tBinOfMaxResponse)-temporalSupportSeconds(midPoint);
    delaySamples = round(theDelaySeconds/ (temporalSupportSeconds(2)-temporalSupportSeconds(1)));
    theCenteredWaveform = circshift(theWaveform, -delaySamples);

    % Window it
    alpha = 0.25;
    win = tukeywin(numel(temporalSupportSeconds), alpha);
    win = reshape(win, size(theCenteredWaveform));

    % Do the windowing
    theCenteredWaveform = theCenteredWaveform .* win;
end