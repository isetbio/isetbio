%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz)

     % gain
    initialValues(1) = 10;
    lowerBounds(1) = 10;
    upperBounds(1) = 200;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 40;
    lowerBounds(numel(lowerBounds)+1) = -50;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.25;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.25;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';


    % lowpass2 gain
    initialValues(numel(initialValues)+1) = 0.9;
    lowerBounds(numel(lowerBounds)+1) = 0.5;
    upperBounds(numel(upperBounds)+1) = 1.5;
    paramNames{numel(paramNames)+1} = 'LP2 filter gain';


    % phase of LP1
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 720;
    paramNames{numel(paramNames)+1} = 'LP1 phase (degs)';

    % phase of LP2
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 720;
    paramNames{numel(paramNames)+1} = 'LP2 phase (degs)';


    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;

    timeConstant1Seconds = theCurrentParams(3)*1e-3;
    timeConstant2Seconds = theCurrentParams(4)*1e-3;
    
    % filter orders must be integer
    keepFilterOrdersInteger = ~true;
    if (keepFilterOrdersInteger)
        theCurrentParams(5) = round(theCurrentParams(5));
        theCurrentParams(6) = round(theCurrentParams(6));
    end
    lowPassFilterOrder = theCurrentParams(5);
    lowPass2FilterOrder = theCurrentParams(6);
    
    lowpass2gain = theCurrentParams(7);


    phaseDegs = mod(theCurrentParams(8),360);
    phaseRadians = phaseDegs/180*pi;

    phaseDegs = mod(theCurrentParams(9),360);
    phaseRadians2 = phaseDegs/180*pi;

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1 + 1i*phaseRadians + 1i * omega * timeConstant1Seconds) .^ (-lowPassFilterOrder);
    theLowPass2FilterTTF = (1 + 1i*phaseRadians2 + 1i * omega * timeConstant2Seconds) .^ (-lowPass2FilterOrder);

    theFilterTTF = gain * theDelayFilterTTF .* (theLowPassFilterTTF - lowpass2gain*theLowPass2FilterTTF); 
end

