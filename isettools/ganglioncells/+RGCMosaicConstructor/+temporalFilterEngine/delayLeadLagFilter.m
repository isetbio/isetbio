%
% RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter
%
%
 
function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz)
    

    % gain
    initialValues(1) = 0.13;
    lowerBounds(1) = 0.05;
    upperBounds(1) = 10;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = -30;
    lowerBounds(numel(lowerBounds)+1) = -60;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % timeConstantSeconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % deltaTau (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 200;
    paramNames{numel(paramNames)+1} = 'HP delta time constant (msec)';
    
    % leadLagFilterOrder
    initialValues(numel(initialValues)+1) = 3;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order';
    

    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.65;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;
    highPassTimeConstantSeconds = theCurrentParams(3)*1e-3;
    deltaTau = theCurrentParams(4)*1e-3;
    
    lowPassTimeConstantSeconds = theCurrentParams(6)*1e-3;
    lowPass2TimeConstantSeconds = theCurrentParams(7)*1e-3;

    % filter orders must be integer
    keepFilterOrdersInteger = false;
    if (keepFilterOrdersInteger)
        theCurrentParams(5) = round(theCurrentParams(5));
        theCurrentParams(8) = round(theCurrentParams(8));
        theCurrentParams(9) = round(theCurrentParams(9));
    end

    leadLagFilterOrder = theCurrentParams(5);
    lowPassFilterOrder = theCurrentParams(8);
    lowPass2FilterOrder = theCurrentParams(9);
    
    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilterTTF = ((1 + 1i * omega * (highPassTimeConstantSeconds+deltaTau)) ./ (1 + 1i * omega * highPassTimeConstantSeconds) ) .^ leadLagFilterOrder;

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-lowPassFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPass2FilterTTF = (1 + 1i * omega * lowPass2TimeConstantSeconds) .^ (-lowPass2FilterOrder);


    theFilterTTF = gain * theDelayFilterTTF .* theLeadLagFilterTTF .* theLowPassFilterTTF .* theLowPass2FilterTTF;
end
