function  velocityTimeSeries = computeVelocity(obj, theEMpath)
% Compute the velocity for the passed emPath
%
% Syntax:
%   velocityTimeSeries = computeVelocity(obj, theEMpath)
%   velocityTimeSeries = obj.computeVelocity(theEMpath)
%
% Description:
%    Compute the velocity for the fixationalEM using the provided emPath
%    using a 3rd order Savitzky-Golay temporal smoothing filter defined 
%    over obj.velocityMeasurementIntervalSeconds. Note that the velocity
%    will be highly dependent on this time window as the emPath is a
%    modified Brownian motion process.
%
% Inputs:
%    obj                - Object. The fixationalEM object.
%    theEMpath          - Matrix. A 2xN matrix containing the movement path
%                         for each eye.
%
% Outputs:
%    velocityTimeSeries - Vector. The velocity time series.
%
% Optional key/value pairs:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
%

if (size(theEMpath, 1) ~= 2), theEMpath = theEMpath'; end
timeStepSeconds = obj.timeAxis(2) - obj.timeAxis(1);
timeSamplesNum = size(theEMpath, 2);
xPos = squeeze(theEMpath(1, :));
yPos = squeeze(theEMpath(2, :));
filterOrder = 3;
filterLengthSamples = round(obj.velocityMeasurementIntervalSeconds / ...
    timeStepSeconds);
if (mod(filterLengthSamples, 2) == 0)
    filterLengthSamples = filterLengthSamples + 1;
end
if (filterOrder<filterLengthSamples)
    xPosFiltered = sgolayfilt(xPos, filterOrder, filterLengthSamples);
    yPosFiltered = sgolayfilt(yPos, filterOrder, filterLengthSamples);
else
    xPosFiltered = xPos;
    yPosFiltered = yPos;
end

velocityTimeSeries = zeros(1, timeSamplesNum);
velocityTimeSeries(2:timeSamplesNum) = sqrt(((diff(xPosFiltered)) / ...
    timeStepSeconds) .^ 2 + ((diff(yPosFiltered)) / timeStepSeconds) .^ 2);

end
