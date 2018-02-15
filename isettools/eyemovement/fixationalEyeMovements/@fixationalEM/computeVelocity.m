function  velocityTimeSeries = computeVelocity(obj, theEMpath)
    if (size(theEMpath,1) ~= 2)
        theEMpath = theEMpath';
    end
    timeStepSeconds = obj.timeAxis(2)-obj.timeAxis(1);
    timeSamplesNum = size(theEMpath,2);
    xPos = squeeze(theEMpath(1,:));
    yPos = squeeze(theEMpath(2,:));
    filterOrder = 3;
    filterLengthSamples = round(obj.velocityMeasurementIntervalSeconds/timeStepSeconds);
    if (mod(filterLengthSamples,2) == 0)
        filterLengthSamples = filterLengthSamples + 1;
    end
    if (filterOrder<filterLengthSamples)
        xPosFiltered = sgolayfilt(xPos, filterOrder, filterLengthSamples);
        yPosFiltered = sgolayfilt(yPos, filterOrder, filterLengthSamples);
    else
        xPosFiltered = xPos;
        yPosFiltered = yPos;
    end
    
    velocityTimeSeries = zeros(1,timeSamplesNum);
    velocityTimeSeries(2:timeSamplesNum) = sqrt(((diff(xPosFiltered))/timeStepSeconds).^2 + ...
                                   ((diff(yPosFiltered))/timeStepSeconds).^2);
end
