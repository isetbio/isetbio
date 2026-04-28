%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz)


    keepFilterOrdersInteger = ~true;

    

    % gain
    initialValues(1) = 100;
    lowerBounds(1) = 1;
    upperBounds(1) = 1000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 25;
    lowerBounds(numel(lowerBounds)+1) = 15;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % lowpass1 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.8;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 10;
    paramNames{numel(paramNames)+1} = 'LP1 time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.9;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 10;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';



    % lowpass1 FilterOrder
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'nLTL-1';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'nLTl-2';

    % differentiation gain
    initialValues(numel(initialValues)+1) = 0.95;
    lowerBounds(numel(lowerBounds)+1) = 0.6;
    upperBounds(numel(upperBounds)+1) = 1.1;
    paramNames{numel(paramNames)+1} = 'differentiation gain';   


    % differentiation delay
    initialValues(numel(initialValues)+1) = 0.;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 0;
    paramNames{numel(paramNames)+1} = 'differentiation delay (msec)';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;

    % Time constants
    timeConstant1Seconds = theCurrentParams(3)*1e-3;
    timeConstant2Seconds = theCurrentParams(4)*1e-3;


    % Filter orders must be integer
    if (keepFilterOrdersInteger)
        theCurrentParams(5) = round(theCurrentParams(5));
        theCurrentParams(6) = round(theCurrentParams(6));
    end

    lowPassFilterOrder = theCurrentParams(5)/theCurrentParams(3);
    lowPassFilter2Order = theCurrentParams(6)/theCurrentParams(4);

    % Differntiation gain
    differentiationGain = theCurrentParams(7);

    % Differntiation delay
    differentiationDelaySeconds = theCurrentParams(8)*1e-3;
   

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1 + 1i * (omega * timeConstant1Seconds)) .^ (-lowPassFilterOrder);
    theLowPassFilter2TTF = (1 + 1i * (omega * timeConstant2Seconds)) .^ (-lowPassFilter2Order);
   
    theFilterTTF = gain * theDelayFilterTTF .* (theLowPassFilterTTF - differentiationGain * exp(-1i * omega * differentiationDelaySeconds) .* theLowPassFilter2TTF);

end

