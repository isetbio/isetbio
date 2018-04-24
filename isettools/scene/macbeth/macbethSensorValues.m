function [sensorImg, sensorSD, cornerPoints] = ...
    macbethSensorValues(sensor, showSelection, cornerPoints)
% Identify MCC patches and calculate RGB mean (s.d.) of the 24 patches
%
% Syntax:
%	[sensorImg, sensorSD, cornerPoints] = ...
%       macbethSensorValues(sensor, showSelection, cornerPoints);
%
% Description:
%    This routine is designed to analyze sensor RGB values. It calculates
%    how to (linearly) transform the sensor RGB into the ideal RGB values
%    for and MCC target. We can find this transform by
%
%       idealLRGB = sensorImg L(3x3)
%
%    The transform can be adjusted to minimize different types of errors,
%    such as the achromatic series, delta E, or RMSE. More routines will be
%    developed for this purpose. If we find the transform L in different
%    circumstances, say different ambient lighting for the MCC, then we can
%    build up a whole set that is appropriate for color balancing.
%
%    There are examples contained within the code. To access, enter 'edit
%    macbethSensorValues.m' into the Command Window.
%
% Inputs:
%    sensor        - (Optional) The sensor. Default is to get a vc Object
%                    of type 'sensor'.
%    showSelection - (Optional) Boolean indicating whether or not to
%                    display the selection. Default is true.
%    cornerPoints  - (Optional) Rectangular Corner points. Following
%                    macbethRectangles, the points are A 4 x 2 matrix of
%                    four points indicating the corners of the region we
%                    want to parcel up. The order is lower left, lower
%                    right, upper right, upper left. They can be chosen
%                    using the GUI and macbethSelect. The entries are (row,
%                    col), which is (y, x).
%
% Outputs:
%    sensorImg     - The sensor image
%    sensorSD      - 
%    cornerPoints  - The same corner points (if provided), else the
%                    calculated corners of the MCC patches
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: showSelection is not working any more.
%    * TODO: Either create macbethSelect, or deprecate this function.
%
% See Also:
%    macbethSelect, macbethCompareIdeal
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/31/18  jnm  Formatting


% Example:
%{
    sensor = vcGetObject('sensor');
    [sensorImg, sensorSD, cornerPoints] = macbethSensorValues(sensor, 1);
%}


if notDefined('sensor'), sensor = vcGetObject('sensor'); end
if notDefined('showSelection'), showSelection = true; end
fullData = true;
% Get the raw sensor data
if notDefined('cornerPoints')
    [fullRGB, ~, ~, cornerPoints] = macbethSelect(sensor, ...
        showSelection, fullData);
else
    [fullRGB, ~, ~, cornerPoints] = macbethSelect(sensor, ...
        showSelection, fullData, cornerPoints);
end

nSensors = size(fullRGB{1}, 2);
sensorImg = zeros(24, nSensors);
if nargout == 2, sensorSD = zeros(24, nSensors); else, sensorSD = []; end

% Fix up the NaNs for the sensor data
for ii = 1:24  % For each chip
    tmp = fullRGB{ii};
    for band = 1:nSensors  % For each band
        foo = tmp(:, band);
        sensorImg(ii, band) = mean(foo(~isnan(foo)));
        if nargout == 2, sensorSD(ii, band) = std(foo(~isnan(foo))); end
    end
end

end