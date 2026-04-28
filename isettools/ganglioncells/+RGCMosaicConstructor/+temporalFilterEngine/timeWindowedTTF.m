%
% RGCMosaicConstructor.temporalFilterEngine.timeWindowedTTF(achievedTTF, temporalWeightingLimitsSeconds)
%

function [theTimeWindowedTTF, theOneSidedTimeWindowedTTF] = timeWindowedTTF(theUncroppedUncenteredTTF, temporalFrequencySupportHz, temporalLimitsSeconds)

    % Generate the corresponding temporal impulse response function 
    theTemporalResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theUncroppedUncenteredTTF, temporalFrequencySupportHz, ...
                'upsample', 1);

    % Find the time of the peak
    [~, tBinOfMaxResponse] = max(abs(theTemporalResponseData.amplitude(:)));
    shiftSamplesToCenterImpulseResponse = 0.5*(numel(theTemporalResponseData.temporalSupportSeconds))-tBinOfMaxResponse;

    % Center the peak of the response
    theTemporalResponseData.amplitude = circshift(theTemporalResponseData.amplitude, shiftSamplesToCenterImpulseResponse);

    % Window it
    dT = theTemporalResponseData.temporalSupportSeconds(2)-theTemporalResponseData.temporalSupportSeconds(1);
    theWindow = hann(round((temporalLimitsSeconds(2)-temporalLimitsSeconds(1))/dT));
    shiftSamples = 0.5*(numel(theTemporalResponseData.temporalSupportSeconds))-0.5*numel(theWindow);
    theFullWindow = theTemporalResponseData.temporalSupportSeconds*0;
    theFullWindow(1:numel(theWindow)) = theWindow;
    theFullWindow = circshift(theFullWindow, shiftSamples);

    theCroppedTemporalImpulseResponseData.amplitude = theFullWindow .* theTemporalResponseData.amplitude;


    % FFT of the centered, cropped IR
    theTimeWindowedTTF = fft(theCroppedTemporalImpulseResponseData.amplitude);


    % Extract the one-sided TTF spectrum of the centered, cropped IR
    % 1. Back to time-domain
    tmp = ifft(theTimeWindowedTTF);
    % 2. Undo the shift
    tmp = circshift(tmp, -shiftSamplesToCenterImpulseResponse);
    % 3. Back to Fourier domain
    tmp = fft(tmp);
    % 4. Take the one sided part
    theOneSidedTimeWindowedTTF = 2*tmp(1:numel(temporalFrequencySupportHz));
    % 5. Reshape to original shape
    theOneSidedTimeWindowedTTF = reshape(theOneSidedTimeWindowedTTF, size(theUncroppedUncenteredTTF));

    
    debug = ~true;
    if (debug)

        theTemporalResponseDataAfterCropping = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    theOneSidedTimeWindowedTTF, temporalFrequencySupportHz);
    
        figure(111); clf
        subplot(1,2,1)
        plot(theTemporalResponseData.temporalSupportSeconds, theTemporalResponseData.amplitude, 'ko');
        hold on;
        plot(theTemporalResponseDataAfterCropping.temporalSupportSeconds, theTemporalResponseDataAfterCropping.amplitude, 'r-', 'LineWidth', 1.0);
        plot(theTemporalResponseData.temporalSupportSeconds, theFullWindow, 'g-');
        subplot(1,2,2);
        plot(temporalFrequencySupportHz, abs(theUncroppedUncenteredTTF), 'ko');
        hold on;
        plot(temporalFrequencySupportHz, abs(theOneSidedTimeWindowedTTF), 'r-', 'LineWidth',1.0);
        set(gca, 'XScale', 'log', 'XLim', [0.3 300], 'XTick', [0.3 1 3 10 30 100 300]);
        pause;
    end
   
end
