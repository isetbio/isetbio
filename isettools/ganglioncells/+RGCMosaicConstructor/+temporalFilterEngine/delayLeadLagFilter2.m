%
% RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter2
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayLeadLagFilter2(theCurrentParams, temporalFrequencySupportHz)
    
    % gain
    initialValues(1) = 0.13;
    lowerBounds(1) = 0.05;
    upperBounds(1) = 10;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = -30;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % timeConstant1Seconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 2.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 25;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % timeConstant2Seconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 2.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 25;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';

    % deltaTau1 (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 300;
    paramNames{numel(paramNames)+1} = 'HP delta time constant 1 (msec)';
    
    % deltaTau2 (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 300;
    paramNames{numel(paramNames)+1} = 'HP delta time constant 2 (msec)';

    % leadLagFilterOrder1
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order 1';
    
    % leadLagFilterOrder2
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order 2';



    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2.0;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';

    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 2;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    
    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;
    highPassTimeConstant1Seconds = theCurrentParams(3)*1e-3;
    highPassTimeConstant2Seconds = theCurrentParams(4)*1e-3;
    deltaTau1 = theCurrentParams(5)*1e-3;
    deltaTau2 = theCurrentParams(6)*1e-3;

    % filter orders must be integer
    theCurrentParams(7) = round(theCurrentParams(7));
    leadLagFilter1Order = theCurrentParams(7);

    % filter orders must be integer
    theCurrentParams(8) = round(theCurrentParams(8));
    leadLagFilter2Order = theCurrentParams(8);


    lowPassTimeConstantSeconds = theCurrentParams(9)*1e-3;
    lowPass2TimeConstantSeconds = theCurrentParams(10)*1e-3;

    
    % filter orders must be integer
    theCurrentParams(11) = round(theCurrentParams(11));
    lowPassFilterOrder = theCurrentParams(11);
   
    % filter orders must be integer
    theCurrentParams(12) = round(theCurrentParams(12));
    lowPass2FilterOrder = theCurrentParams(12);
    
    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilter1TTF = ((1 + 1i * omega * (highPassTimeConstant1Seconds+deltaTau1)) ./ (1 + 1i * omega * highPassTimeConstant1Seconds) ) .^ round(leadLagFilter1Order);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilter2TTF = ((1 + 1i * omega * (highPassTimeConstant2Seconds+deltaTau2)) ./ (1 + 1i * omega * highPassTimeConstant2Seconds) ) .^ round(leadLagFilter2Order);


    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-lowPassFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPass2FilterTTF = (1 + 1i * omega * lowPass2TimeConstantSeconds) .^ (-lowPass2FilterOrder);


    theFilterTTF = gain * theDelayFilterTTF .* theLeadLagFilter1TTF .* theLeadLagFilter2TTF .* theLowPassFilterTTF .* theLowPass2FilterTTF;
end

