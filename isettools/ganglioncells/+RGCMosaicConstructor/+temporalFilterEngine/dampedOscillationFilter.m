%
% RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter
%
%

%{
%   Some good parameters
    gain: 1400.58
    delay (msec): 15
    a1 (oscillation frequency): 115.607
    a2 (modulation of frequency over time): 37.1896
    a3 (steepness of exponential decay): 222.208
    phase (phase of oscillation): -5.82764


    theFrequencySupportHz = 0:0.5:200;
    [~, initialValues] = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter([],[]);
    theFilterTTF = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(...
        initialValues, theFrequencySupportHz);

    theImpulseResponseFunctionStruct = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theFilterTTF, theFrequencySupportHz, ...
        'causal', false);

    idx = find(theImpulseResponseFunctionStruct.temporalSupportSeconds<0.25);
    theImpulseResponseFunctionStruct.temporalSupportSeconds = theImpulseResponseFunctionStruct.temporalSupportSeconds(idx);
    theImpulseResponseFunctionStruct.amplitude = theImpulseResponseFunctionStruct.amplitude(idx);

    theFrequencySupportHz = 0.5:0.5:200;
    nWarmpUpCycles = 4;
    dtSeconds = theImpulseResponseFunctionStruct.temporalSupportSeconds(2)-theImpulseResponseFunctionStruct.temporalSupportSeconds(1)
    for iTF = 1:numel(theFrequencySupportHz)
        temporalFrequencyHz = theFrequencySupportHz(iTF);
        tOneStimulusCycle = 1/temporalFrequencyHz
        temporalSupportSeconds = 0:dtSeconds:(nWarmpUpCycles*tOneStimulusCycle);
        theSinusoidalInput = sin(2*pi*temporalFrequencyHz*temporalSupportSeconds);
      
        [theSinusoidalInput(1) theSinusoidalInput(end)]
        theFullResponse = conv(theSinusoidalInput(:),theImpulseResponseFunctionStruct.amplitude(:));
        theFullResponse = theFullResponse(1:numel(temporalSupportSeconds));
        temporalSupportSecondsFullResponse = temporalSupportSeconds;

        theFullResponse = 0.5*theFullResponse / max(abs(theFullResponse));

        
        figure(1); clf;
        subplot(3,1,1);
        % plot the impulse response
        plot(theImpulseResponseFunctionStruct.temporalSupportSeconds, theImpulseResponseFunctionStruct.amplitude, 'r-');

        
        % Plot the stimulus and the full response
        subplot(3,1,2);
        plot(temporalSupportSeconds, theSinusoidalInput, 'k-');
        hold on
        plot(temporalSupportSecondsFullResponse, theFullResponse, 'r-');

        % Trim the stimulus response to the last period
        idx = find(temporalSupportSecondsFullResponse >= temporalSupportSecondsFullResponse(end)-tOneStimulusCycle);
        temporalSupportSecondsLastPeriod = temporalSupportSecondsFullResponse(idx);

        t1 = temporalSupportSecondsLastPeriod(1);
        temporalSupportSecondsLastPeriod = temporalSupportSecondsLastPeriod-t1;
        theResponse = theFullResponse(idx);

        idxSinusoidalInput = find(temporalSupportSeconds>=t1);
        theStimulus = theSinusoidalInput(idxSinusoidalInput);
        temporalSupportSeconds = temporalSupportSeconds(idxSinusoidalInput);
        temporalSupportSeconds = temporalSupportSeconds-t1;

        subplot(3, 1,3);
        plot(temporalSupportSecondsLastPeriod,theResponse, 'r-');
        hold on;
        plot(temporalSupportSeconds,theStimulus, 'k-');
        
        pause
    end




    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1000 1000]);
    ax = subplot(2,1,1);
    plot(ax,theImpulseResponseFunctionStruct.temporalSupportSeconds*1e3, ...
        theImpulseResponseFunctionStruct.amplitude, ...
        'r-', 'LineWidth', 1.5);
    set(ax, 'XLim', [0 300]);
    ylabel(ax, 'response')
    xlabel(ax, 'time (msec)')
    set(ax, 'FontSize', 16)
    grid(ax, 'on')

    ax = subplot(2,1,2);
    yyaxis(ax, 'left')
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, theFrequencySupportHz, ...
        theFilterTTF, '-', ...
        true, false, [1 0 0], ...
        '');

    yyaxis(ax, 'right')
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, theFrequencySupportHz, ...
        theFilterTTF, '-', ...
        true, false, [0 0 1], ...
        '');
    grid(ax, 'on')

%}


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    dampedOscillationFilter(theCurrentParams, temporalFrequencySupportHz)

    % gain
    initialValues(1) = 1.0;
    lowerBounds(1) = 0.01;
    upperBounds(1) = 10.0;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % a1
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 120;
    paramNames{numel(paramNames)+1} = 'a1 (oscillation frequency)';

    % a2
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a2 (modulation of frequency over time)';

    % a3
    initialValues(numel(initialValues)+1) = 130;
    lowerBounds(numel(lowerBounds)+1) = 50;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a3 (steepness of exponential decay)';

    % phase of oscillation
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = -180;
    upperBounds(numel(upperBounds)+1) = 180;
    paramNames{numel(paramNames)+1} = 'phase (phase of oscillation)';



    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;
    a1 = theCurrentParams(3);
    a2 = theCurrentParams(4);
    a3 = theCurrentParams(5);
    phaseDegs = theCurrentParams(6);

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);


    M = numel(temporalFrequencySupportHz);
    N  = 2 * (M - 1);
    samplingFrequency = 2 * temporalFrequencySupportHz(end);
    temporalSupportSeconds = (0:N-1) / samplingFrequency;


    tCausal = max(0,temporalSupportSeconds-delaySeconds);

    dampedOscillationImpulseResponse = temporalSupportSeconds .* ...
        sin(phaseDegs/180*pi + 2*pi*(a1*tCausal .* (tCausal+1).^(-a2))) .* ...
        exp(-a3*tCausal);


    dampedOscillationTTFtwoSided = fft(dampedOscillationImpulseResponse);

    % Generate one sided TTF
    dampedOscillationTTFoneSided(1:M) = dampedOscillationTTFtwoSided(1:M);
    dampedOscillationTTFoneSided(2:M) = 2*dampedOscillationTTFtwoSided(2:M);

    theFilterTTF = 1e3 * gain * theDelayFilterTTF .* dampedOscillationTTFoneSided;
end
