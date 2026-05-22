%
% RGCMosaicConstructor.temporalFilterEngine.asymmetricBandPassFilter
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    asymmetricBandPassFilter(theCurrentParams, temporalFrequencySupportHz)

    % gain
    initialValues(1) = 8000;
    lowerBounds(1) = 1;
    upperBounds(1) = 20000;
    paramNames{1} = 'gain';

    % delaySeconds
    initialValues(numel(initialValues)+1) = 40;
    lowerBounds(numel(lowerBounds)+1) = -10;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % tau lowpass
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 3;
    paramNames{numel(paramNames)+1} = 'tau lowpass (msec)';


    % tau Low / tau High difference (must be > 1)
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'tau lowpass / tau highpass(msec)';


    % order (lowpass) (must be >= 0)
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'order lowpass';

    % order (highpass) (must be >= 1)
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'order highpass';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end

 
    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;

    tauLowPassSeconds = theCurrentParams(3)*1e-3;
    tauHighPassSeconds = max([1 theCurrentParams(3)*theCurrentParams(4)])*1e-3;
    
    % WE MUST HAVE: tauHighPassSeconds < tauLowPassSeconds,
    % so that highPassFrequency > lowPassFrequency

    keepFilterOrdersInteger = ~true;
    if (keepFilterOrdersInteger)
        theCurrentParams(5) = max([0 round(theCurrentParams(5))]);
        theCurrentParams(6) = max([1 round(theCurrentParams(6))]);
    end

    orderLowPass = theCurrentParams(5);
    orderHighPass = theCurrentParams(6);

    highPassFrequency = 1/tauHighPassSeconds;
    lowPassFrequency = 1/tauLowPassSeconds;
    centerFrequency = sqrt(highPassFrequency * lowPassFrequency);


    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1/tauLowPassSeconds ./ (1i * omega + 1/tauLowPassSeconds)) .^ orderLowPass;
    theHighPassFilterTTF = ((1i * omega) ./ (1i * omega + 1/tauHighPassSeconds)) .^ orderHighPass;
   

    theFilterTTF = gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;

end
