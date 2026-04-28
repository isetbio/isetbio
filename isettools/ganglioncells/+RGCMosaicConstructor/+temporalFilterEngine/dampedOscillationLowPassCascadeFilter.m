%
% RGCMosaicConstructor.temporalFilterEngine.dampedOscillationLowPassCascadeFilter
%
%

%{
%   Some good parameters
    gain: 14.9942
    delay (msec): 25
    a1 (oscillation frequency): 110.683
    a2 (modulation of frequency over time): 23.206
    a3 (steepness of exponential decay): 272.017
    phase (phase of oscillation): 203.047
    LP time constant (msec): 0.750553
    nStages x Tau (LP): 28.5834


    theFrequencySupportHz = 0:0.5:300;
    [~, initialValues] = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter([],[]);
    theFilterTTF = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(...
        initialValues, theFrequencySupportHz);

    theImpulseResponseFunctionStruct = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theFilterTTF, theFrequencySupportHz, ...
        'causal', false);

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
    dampedOscillationLowPassCascadeFilter(theCurrentParams, temporalFrequencySupportHz)

    % gain
    initialValues(1) = 4.0;
    lowerBounds(1) = 1.0;
    upperBounds(1) = 50.0;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 25;
    lowerBounds(numel(lowerBounds)+1) = 25;
    upperBounds(numel(upperBounds)+1) = 25;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % a1
    initialValues(numel(initialValues)+1) = 80;
    lowerBounds(numel(lowerBounds)+1) = 50;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'a1 (oscillation frequency)';

    % a2
    initialValues(numel(initialValues)+1) = 15;
    lowerBounds(numel(lowerBounds)+1) = 5;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'a2 (modulation of frequency over time)';

    % a3
    initialValues(numel(initialValues)+1) = 215;
    lowerBounds(numel(lowerBounds)+1) = 100;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a3 (steepness of exponential decay)';

    % phase of oscillation
    initialValues(numel(initialValues)+1) = 210;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'phase (phase of oscillation)';


    % lowpas timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.36;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';


    % Product of filter orded x time constant (LP)
    initialValues(numel(initialValues)+1) = 25;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'nStages x Tau (LP)';


    lowerBounds = initialValues*(1-0.5);
    upperBounds = initialValues*(1+0.5);

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

    lowPassFilterTimeConstantSeconds = theCurrentParams(7)*1e-3;
    lowPassFilterOrder = theCurrentParams(8)/theCurrentParams(7);



    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);


    M = numel(temporalFrequencySupportHz);
    N  = 2 * (M - 1);
    samplingFrequency = 2 * temporalFrequencySupportHz(end);
    temporalSupportSeconds = (0:N-1) / samplingFrequency;


    dampedOscillationImpulseResponse = temporalSupportSeconds .* ...
        sin(phaseDegs/180*pi + 2*pi*(a1*temporalSupportSeconds .* (abs(temporalSupportSeconds)+1).^(-a2))) .* ...
        exp(-a3*temporalSupportSeconds);

    % Center the waveform
    [theCenteredIR, theCenteringDelaySeconds] = ...
        RGCMosaicConstructor.temporalFilterEngine.centerAndWindowTemporalImpulseResponse(temporalSupportSeconds, dampedOscillationImpulseResponse);

    dampedOscillationTTFtwoSided = fft(theCenteredIR);

    % Generate one sided TTF
    dampedOscillationTTFoneSided(1:M) = dampedOscillationTTFtwoSided(1:M);
    dampedOscillationTTFoneSided(2:M) = 2*dampedOscillationTTFtwoSided(2:M);

    % Undo the centering
    dampedOscillationTTFoneSided = exp(-1i * omega * (-theCenteringDelaySeconds)) .* dampedOscillationTTFoneSided;

    % Synthesize the low-pass filter
    theLowPassFilterTTF = (1 + 1i * (omega * lowPassFilterTimeConstantSeconds)) .^ (-lowPassFilterOrder);

    theFilterTTF = 1e4 * gain * theDelayFilterTTF .* theLowPassFilterTTF .* dampedOscillationTTFoneSided;
    
end



