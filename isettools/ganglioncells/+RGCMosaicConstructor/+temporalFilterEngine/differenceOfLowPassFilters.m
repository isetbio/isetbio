%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters
%
%

%{

gain: 999.989
delay (msec): 15
LP1 time constant (msec): 0.378987
LP2 time constant (msec): 0.614186
nLTL-1: 36.8432
nLTl-2: 36.5348
differentiation gain: 1.00659
differentiation delay (msec): 0

%}


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz)


    keepFilterOrdersInteger = ~true;

    
    % gain
    initialValues(1) = 1000;
    lowerBounds(1) = 500;
    upperBounds(1) = 2000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 15;
    lowerBounds(numel(lowerBounds)+1) = 15;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % lowpass1 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.37;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'LP1 time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.61;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';



    % Product of filter orded x time constant (LP-1)
    initialValues(numel(initialValues)+1) = 36.8;
    lowerBounds(numel(lowerBounds)+1) = 20;
    upperBounds(numel(upperBounds)+1) = 45;
    paramNames{numel(paramNames)+1} = 'nStages x Tau (LP-1)';

     % Product of filter orded x time constant (LP-2)
    initialValues(numel(initialValues)+1) = 36.5;
    lowerBounds(numel(lowerBounds)+1) = 20;
    upperBounds(numel(upperBounds)+1) = 45;
    paramNames{numel(paramNames)+1} = 'nStages x Tau (LP-1)';

    % differentiation gain
    initialValues(numel(initialValues)+1) = 0.95;
    lowerBounds(numel(lowerBounds)+1) = 0.6;
    upperBounds(numel(upperBounds)+1) = 1.0;
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

