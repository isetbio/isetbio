%
% RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters
%
%

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz)


    cascadeWithThirdLowPassFilter = false;


     % gain
    initialValues(1) = 60;
    lowerBounds(1) = 1;
    upperBounds(1) = 200;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 150;
    paramNames{numel(paramNames)+1} = 'delay (msec)';


    % lowpass1 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 15;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP1 time constant (msec)';


    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 12;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';



    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'LP1 filter order';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 40;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';


 
    % phase of LP1
    initialValues(numel(initialValues)+1) = 16;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'LP1 phase (degs)';

    % phase of LP2
    initialValues(numel(initialValues)+1) = 31;
    lowerBounds(numel(lowerBounds)+1) = -360;
    upperBounds(numel(upperBounds)+1) = 360;
    paramNames{numel(paramNames)+1} = 'LP2 phase (degs)';


    % differentiation gain
    initialValues(numel(initialValues)+1) = 0.4;
    lowerBounds(numel(lowerBounds)+1) = 0.01;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'differentiation gain';


   
    
    if (cascadeWithThirdLowPassFilter)
        % lowpass3 timeConstantSeconds
        initialValues(numel(initialValues)+1) = 5;
        lowerBounds(numel(lowerBounds)+1) = 0.01;
        upperBounds(numel(upperBounds)+1) = 20;
        paramNames{numel(paramNames)+1} = 'LP3 time constant (msec)';
        
        % lowpass3 FilterOrder
        initialValues(numel(initialValues)+1) = 0; %11;
        lowerBounds(numel(lowerBounds)+1) = 0; %1;
        upperBounds(numel(upperBounds)+1) = 0; %20;
        paramNames{numel(paramNames)+1} = 'LP3 filter order';
    
    
        % phase of LP3
        initialValues(numel(initialValues)+1) = 31;
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
    keepFilterOrdersInteger = ~true;
    if (keepFilterOrdersInteger)
        theCurrentParams(5) = round(theCurrentParams(5));
        theCurrentParams(6) = round(theCurrentParams(6));
        
    end
    lowPassFilterOrder = theCurrentParams(5);
    lowPass2FilterOrder = theCurrentParams(6);


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
            theCurrentParams(8) = round(theCurrentParams(11));
            lowPass3FilterOrder = theCurrentParams(11);
        end

         phaseDegs = theCurrentParams(12);
         phaseRadians3 = phaseDegs/180*pi;
    end
   

    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1 + 1i * (omega * timeConstant1Seconds - phaseRadians)) .^ (-lowPassFilterOrder);
    theLowPass2FilterTTF = (1 + 1i * (omega * timeConstant2Seconds - phaseRadians2)) .^ (-lowPass2FilterOrder);
   

    theFilterTTF = gain * theDelayFilterTTF .* (theLowPassFilterTTF -  differentiationGain * theLowPass2FilterTTF);
    
    if (cascadeWithThirdLowPassFilter)
        theLowPass3FilterTTF = (1 + 1i * (omega * timeConstant3Seconds - phaseRadians3)) .^ (-lowPass3FilterOrder);
        theFilterTTF = theFilterTTF.* theLowPass3FilterTTF; 
    end

end

