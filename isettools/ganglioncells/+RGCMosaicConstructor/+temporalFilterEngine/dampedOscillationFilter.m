%
% RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    dampedOscillationFilter(theCurrentParams, temporalFrequencySupportHz)

     % gain
    initialValues(1) = 8000;
    lowerBounds(1) = 1;
    upperBounds(1) = 20000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 32;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 150;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % a1
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 20;
    upperBounds(numel(upperBounds)+1) = 120;
    paramNames{numel(paramNames)+1} = 'a1 (oscillation frequency)';

    % a2
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'a2 (modulation of frequency over time)';

    % a3
    initialValues(numel(initialValues)+1) = 150;
    lowerBounds(numel(lowerBounds)+1) = 50;
    upperBounds(numel(upperBounds)+1) = 400;
    paramNames{numel(paramNames)+1} = 'a3 (steepness of exponential decay)';

    % phase of oscillation
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 0;
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
