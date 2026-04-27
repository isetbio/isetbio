%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz)


    cascadeWithThirdLowPassFilter = ~true;
    keepFilterOrdersInteger = ~true;

    % gain
    initialValues(1) = 220;
    lowerBounds(1) = 1;
    upperBounds(1) = 2000;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 20;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % lowpass1 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 3;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 5;
    paramNames{numel(paramNames)+1} = 'LP1 time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 3.1;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 10;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';



    % lowpass1 FilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 3;
    upperBounds(numel(upperBounds)+1) = 10;
    paramNames{numel(paramNames)+1} = 'LP1 filter order';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 7;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';


 
    % phase of LP1
    initialValues(numel(initialValues)+1) = 40;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'LP1 phase (degs)';

    % phase of LP2
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'LP2 phase (degs)';


    % differentiation gain
    initialValues(numel(initialValues)+1) = 0.9;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 1.5;
    paramNames{numel(paramNames)+1} = 'differentiation gain';


   
    
    if (cascadeWithThirdLowPassFilter)
        % lowpass3 timeConstantSeconds
        initialValues(numel(initialValues)+1) = 12;
        lowerBounds(numel(lowerBounds)+1) = 5;
        upperBounds(numel(upperBounds)+1) = 20;
        paramNames{numel(paramNames)+1} = 'LP3 time constant (msec)';
        
        % lowpass3 FilterOrder
        initialValues(numel(initialValues)+1) = 1;
        lowerBounds(numel(lowerBounds)+1) = 1;
        upperBounds(numel(upperBounds)+1) = 2;
        paramNames{numel(paramNames)+1} = 'LP3 filter order';
    
    
        % phase of LP3
        initialValues(numel(initialValues)+1) = 300;
        lowerBounds(numel(lowerBounds)+1) = -360;
        upperBounds(numel(upperBounds)+1) = 360;
        paramNames{numel(paramNames)+1} = 'LP3 phase (degs)';
    end



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
    lowPassFilterOrder = theCurrentParams(5);
    lowPassFilter2Order = theCurrentParams(6);


    % Phases
    phaseDegs = theCurrentParams(7);
    phaseRadians = phaseDegs/180*pi;

    phaseDegs = theCurrentParams(8);
    phaseRadians2 = phaseDegs/180*pi;

    % Differntiation gain
    differentiationGain = theCurrentParams(9);

    if (cascadeWithThirdLowPassFilter)

        timeConstant3Seconds = theCurrentParams(10)*1e-3;

        if (keepFilterOrdersInteger)
            theCurrentParams(11) = round(theCurrentParams(11));
        end
        lowPassFilter3Order = theCurrentParams(11);

        phaseDegs = theCurrentParams(12);
        phaseRadians3 = phaseDegs/180*pi;
    end
   

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1 + 1i * (omega * timeConstant1Seconds - phaseRadians)) .^ (-lowPassFilterOrder);
    theLowPassFilter2TTF = (1 + 1i * (omega * timeConstant2Seconds - phaseRadians2)) .^ (-lowPassFilter2Order);
   

    theFilterTTF = gain * theDelayFilterTTF .* (theLowPassFilterTTF - differentiationGain * theLowPassFilter2TTF);
    
    if (cascadeWithThirdLowPassFilter)
        theLowPassFilter3TTF = (1 + 1i * (omega * timeConstant3Seconds - phaseRadians3)) .^ (-lowPassFilter3Order);
        theFilterTTF = theFilterTTF.* theLowPassFilter3TTF; 
    end

end

