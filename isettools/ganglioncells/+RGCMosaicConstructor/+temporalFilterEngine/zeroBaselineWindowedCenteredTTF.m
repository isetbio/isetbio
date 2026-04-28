%
% RGCMosaicConstructor.temporalFilterEngine.zeroBaselineWindowedCenteredTTF
%

function [theComplexTTF, delaySeconds, theComplexTTFnonCentered] = ...
    zeroBaselineWindowedCenteredTTF(theComplexTTF, temporalFrequencySupportHz, ...
    minimumDelaySecondsForEstimationOfBaseline, centerIR, verifyOffsetCorrection)

    % Save the original shape, so we can replicate it
    originalShape = size(theComplexTTF);

    % Compute the corresponding impulse response from thePhotocurrentBasedMRGCcellTTF 
    % Here centered in the middle of the time window
    theImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
         theComplexTTF, temporalFrequencySupportHz, ...
         'upsample', 1, ...
         'causal', true);


    % Estimate the baseline from time bins > 0.8 minimumDelaySecondsForEstimationOfBaselin
    idx = find(abs(theImpulseResponseData.temporalSupportSeconds)>=minimumDelaySecondsForEstimationOfBaseline);
    if (isempty(idx))
        error('There are no time samples in the IR (%d  msec) to estimate its baserate for delays < %d msec', ...
            theImpulseResponseData.temporalSupportSeconds(end)*1e3, minimumDelaySecondsForEstimationOfBaseline*1e3);
    end
    baselineOffset = mean(theImpulseResponseData.amplitude(idx));


    % Compute the corresponding impulse response 
    % Not centered, the original one
    theImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
         theComplexTTF, temporalFrequencySupportHz, ...
         'upsample', 1, ...
         'causal', false);

    % Remove that offset
    theImpulseResponseData.amplitude = theImpulseResponseData.amplitude - baselineOffset;

    if (centerIR)
        % Find the time of its peak
        [~, tBinOfMaxResponse] = max(abs(theImpulseResponseData.amplitude(:)));
    
        % Center the impulse response on the time window
        midPoint = round(numel(theImpulseResponseData.temporalSupportSeconds)/2);
        delaySeconds = theImpulseResponseData.temporalSupportSeconds(tBinOfMaxResponse)-theImpulseResponseData.temporalSupportSeconds(midPoint);
        delaySamples = round(delaySeconds / (theImpulseResponseData.temporalSupportSeconds(2)-theImpulseResponseData.temporalSupportSeconds(1)));
        theImpulseResponseData.amplitude = circshift(theImpulseResponseData.amplitude, -delaySamples);
    else
        delaySeconds = 0;
    end


    % Window the impulse response uwing a Tuckey window
    % A Tukey window is very flat in the time domain, close to a value of 1.0 for a majority of the window. 
    % A ‘Taper Length’ can be specified which determines the amount of time that the window maintains 
    % a value of one. The lower the taper length, the longer the Tukey has a value of one over the measurement time.
    alpha = 0.25;
    win = tukeywin(numel(theImpulseResponseData.temporalSupportSeconds), alpha);

    % Do the windowing
    theImpulseResponseData.amplitude = theImpulseResponseData.amplitude .* win;

    % Back to Fourier Domain
    theComplexTTF = fft(theImpulseResponseData.amplitude);

    % Back to one-sided Fourier Domain
    theComplexTTF = theComplexTTF(1:numel(temporalFrequencySupportHz));
    theComplexTTF(2:end) = 2*theComplexTTF(2:end);

    % Back to original shape
    theComplexTTF = reshape(theComplexTTF, originalShape);


    % Undo the delays we introduced to center the two TTFs
    delaySeconds = -delaySeconds;
    omega = 2 * pi * temporalFrequencySupportHz;
    theComplexTTFnonCentered = exp(-1i * omega * delaySeconds) .* theComplexTTF;
    
    if (verifyOffsetCorrection)

        theImpulseResponseDataCentered = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    theComplexTTF, temporalFrequencySupportHz);

        theImpulseResponseDataNonCentered = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    theComplexTTFnonCentered, temporalFrequencySupportHz);
        figure(1111); clf;
        ax = subplot(1,2,1);
        p1 = plot(theImpulseResponseDataCentered.temporalSupportSeconds*1e3, theImpulseResponseDataCentered.amplitude, 'r-');
         xlabel('time (msec)')
        set(gca, 'FontSize', 20)
        ax = subplot(1,2,2);
        p2 = plot(theImpulseResponseDataNonCentered.temporalSupportSeconds*1e3, theImpulseResponseDataNonCentered.amplitude, 'b-');
        xlabel('time (msec)')
        legend([p1 p2], {'centered (fist return arg)', 'non-centered (last return arg)'})
        set(gca, 'FontSize', 20)
        pause
    end

end
