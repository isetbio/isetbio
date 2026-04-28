%
% RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter
%
%

%{
% Some good parameters
gain: 83129.2
delay (msec): 39.8245
a1 (oscillation frequency): 300
a2 (modulation of frequency over time): 99.4001
a3 (steepness of exponential decay): 227.96
phase (phase pof oscillation): -356.539



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
    dampedOscillationFilter(theCurrentParams, temporalFrequencySupportHz)

    % gain
    initialValues(1) = 50000;
    lowerBounds(1) = 1000;
    upperBounds(1) = 100000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 5;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % a1
    initialValues(numel(initialValues)+1) = 120;
    lowerBounds(numel(lowerBounds)+1) = 40;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a1 (oscillation frequency)';

    % a2
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 15;
    upperBounds(numel(upperBounds)+1) = 200;
    paramNames{numel(paramNames)+1} = 'a2 (modulation of frequency over time)';

    % a3
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 50;
    upperBounds(numel(upperBounds)+1) = 500;
    paramNames{numel(paramNames)+1} = 'a3 (steepness of exponential decay)';

    % phase of oscillation
    initialValues(numel(initialValues)+1) = -30;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'phase (phase pof oscillation)';


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


    dampedOscillationInTime = temporalSupportSeconds .* ...
        sin(phaseDegs/180*pi + 2*pi*(a1*temporalSupportSeconds .* (temporalSupportSeconds+1).^(-a2))) .* ...
        exp(-a3*temporalSupportSeconds);

    dampedOscillationTTF = fft(dampedOscillationInTime);
    dampedOscillationTTF = dampedOscillationTTF(1:M);

    theFilterTTF = gain * theDelayFilterTTF .* dampedOscillationTTF;
    
end
