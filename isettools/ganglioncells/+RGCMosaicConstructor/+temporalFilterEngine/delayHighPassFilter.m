%
% RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilter
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayHighPassFilter(theCurrentParams, temporalFrequencySupportHz)
    
    % gain
    initialValues(1) = 100;
    lowerBounds(1) = 0.1;
    upperBounds(1) = 1000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % highpass gain
    initialValues(numel(initialValues)+1) = 0.4;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 10.0;
    paramNames{numel(paramNames)+1} = 'highpaass gain';
  
    % highpass timeConstant msec
    initialValues(numel(initialValues)+1) = 0.01;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 10;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    

    % lowpass time cosntat msec
    initialValues(numel(initialValues)+1) = 0.1;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';
    
    % nL-tL
    initialValues(numel(initialValues)+1) = 40;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'nL-TL (LP)';
    

    % high pass FilterOrder
    initialValues(numel(initialValues)+1) = 0.1;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'nL-TL (HP)';



    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    % Get params
    gain = theCurrentParams(1);                              % Responsitivity at 0 Hz
    conductionDelaySeconds = theCurrentParams(2)*1e-3;       % D
    highPassGain = theCurrentParams(3);                      % Hs
    highPassTimeConstantSeconds = theCurrentParams(4)*1e-3;  % Tau_s (Benardete and Kaplan (1992a) varied this for different contrast levels)
    lowPassTimeConstantSeconds = theCurrentParams(5)*1e-3;   % Tau_l
    
   

    nLowPassStagesNum = theCurrentParams(6)/theCurrentParams(5);               % NlTl
    nHighPassStagesNum = theCurrentParams(7)/theCurrentParams(4);                % always 1 in Benardete & Kaplan (1992a)


    % Circular frequency in radians
    omega = 2 * pi * temporalFrequencySupportHz;

    % Delay filter
    theDelayFilterTTF = exp(-1i * omega * conductionDelaySeconds);

    % 1-stage high-pass filter
    theHighPassFilterTTF = (1 - highPassGain * (1 + 1i * omega * highPassTimeConstantSeconds) .^ (-nHighPassStagesNum));

    % N-stage low-pass filter
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-nLowPassStagesNum);

    theFilterTTF = gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;


end

