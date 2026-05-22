%
% RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter
%
%
 
%{

gain = 21.52
delay = 
HP tau: 3.46
HP delta tau: 4.31
lead-lag filter order: 5.00
LP-1 tau: 0.486
nL-tL-1: 28.33
nL-tL-2: 0.00
%}


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz)
    


    
    % gain
    initialValues(1) = 0.05;
    lowerBounds(1) = 0.01;
    upperBounds(1) = 2;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 15;
    lowerBounds(numel(lowerBounds)+1) = 25;
    upperBounds(numel(upperBounds)+1) = 45;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % timeConstantSeconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 0.11;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % deltaTau (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 10.31;
    lowerBounds(numel(lowerBounds)+1) = 0.3;
    upperBounds(numel(upperBounds)+1) = 15.5;
    paramNames{numel(paramNames)+1} = 'HP delta time constant (msec)';
    
    % leadLagFilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 17;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order';
    

    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.3;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 0.5;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.5;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 1.0;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 28.33;
    lowerBounds(numel(lowerBounds)+1) = 27;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'nLTL LP1';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'nLTL LP1';

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


    leadLagFilterOrder = theCurrentParams(5);

    lowPassFilterOrder = theCurrentParams(8)/theCurrentParams(6);
    lowPass2FilterOrder = theCurrentParams(9)/theCurrentParams(7);
    
    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilterTTF = ((1 + 1i * omega * (highPassTimeConstantSeconds+deltaTau)) ./ (1 + 1i * omega * highPassTimeConstantSeconds)) .^ leadLagFilterOrder;

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-lowPassFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPass2FilterTTF = (1 + 1i * omega * lowPass2TimeConstantSeconds) .^ (-lowPass2FilterOrder);


    theFilterTTF = gain * theDelayFilterTTF .* theLeadLagFilterTTF .* theLowPassFilterTTF .* theLowPass2FilterTTF;
end
