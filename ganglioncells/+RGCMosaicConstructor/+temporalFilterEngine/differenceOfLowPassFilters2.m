%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters2
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters2(theCurrentParams, temporalFrequencySupportHz)

    

    % gain
    initialValues(1) = 100;
    lowerBounds(1) = 1;
    upperBounds(1) = 1000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 35;
    lowerBounds(numel(lowerBounds)+1) = 35;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % low-frequecy (Hz) -3dB
    initialValues(numel(initialValues)+1) = 0.5*200;
    lowerBounds(numel(lowerBounds)+1) = 0.5*200;
    upperBounds(numel(upperBounds)+1) = 0.5*200;
    paramNames{numel(paramNames)+1} = 'low TF (Hz)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 1*200;
    lowerBounds(numel(lowerBounds)+1) = 1*200;
    upperBounds(numel(upperBounds)+1) = 1*200;
    paramNames{numel(paramNames)+1} = 'high TF (Hz)';



    % filterOrder
    initialValues(numel(initialValues)+1) = 4*2;
    lowerBounds(numel(lowerBounds)+1) = 4*2;
    upperBounds(numel(upperBounds)+1) = 4*2;
    paramNames{numel(paramNames)+1} = 'fiter order';

 

    % differentiation gain
    initialValues(numel(initialValues)+1) = 1.0;
    lowerBounds(numel(lowerBounds)+1) = 1.0;
    upperBounds(numel(upperBounds)+1) = 1.0;
    paramNames{numel(paramNames)+1} = 'differentiation gain';   


    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;

    % Low & Hi TF
    f_lowHz = theCurrentParams(3);
    f_highHz = theCurrentParams(4);

 
    filterOrder = theCurrentParams(5);
    

    % Differntiation gain
    differentiationGain = theCurrentParams(6);

    w_low = 2*pi*f_lowHz;
    w_high = 2*pi*f_highHz;

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterLowTTF = (1 + 1i * omega/w_high) .^ (-filterOrder);
    theLowPassFilterHighTTF = (1 + 1i * omega/w_low) .^ (-filterOrder);
   

    theFilterTTF = gain * theDelayFilterTTF .* (theLowPassFilterLowTTF - differentiationGain  * theLowPassFilterHighTTF);

end

