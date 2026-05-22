%
% RGCMosaicConstructor.temporalFilterEngine.dampedOscillationLowPassCascadeFilter
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
    LP time constant (msec): 0.617673
    nStages x Tau (LP): 19.153


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
    initialValues(1) = 10;
    lowerBounds(1) = 1;
    upperBounds(1) = 2000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % a1
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 30;
    upperBounds(numel(upperBounds)+1) = 120;
    paramNames{numel(paramNames)+1} = 'a1 (oscillation frequency)';

    % a2
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 200;
    paramNames{numel(paramNames)+1} = 'a2 (modulation of frequency over time)';

    % a3
    initialValues(numel(initialValues)+1) = 130;
    lowerBounds(numel(lowerBounds)+1) = 100;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a3 (steepness of exponential decay)';

    % phase of oscillation
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = -180;
    upperBounds(numel(upperBounds)+1) = 180;
    paramNames{numel(paramNames)+1} = 'phase (phase of oscillation)';


    % lowpas timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.1;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';


    % Product of filter orded x time constant (LP)
    initialValues(numel(initialValues)+1) = 25;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 80;
    paramNames{numel(paramNames)+1} = 'nStages x Tau (LP)';



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


    tCausal = max(0,temporalSupportSeconds-delaySeconds);

    dampedOscillationImpulseResponse = temporalSupportSeconds .* ...
        sin(phaseDegs/180*pi + 2*pi*(a1*tCausal .* (tCausal+1).^(-a2))) .* ...
        exp(-a3*tCausal);



    dampedOscillationTTFtwoSided = fft(dampedOscillationImpulseResponse);

    % Generate one sided TTF
    dampedOscillationTTFoneSided(1:M) = dampedOscillationTTFtwoSided(1:M);
    dampedOscillationTTFoneSided(2:M) = 2*dampedOscillationTTFtwoSided(2:M);


    % Synthesize the low-pass filter
    theLowPassFilterTTF = (1 + 1i * (omega * lowPassFilterTimeConstantSeconds)) .^ (-lowPassFilterOrder);

    theFilterTTF = gain * theDelayFilterTTF .* theLowPassFilterTTF .* dampedOscillationTTFoneSided;
    
end



