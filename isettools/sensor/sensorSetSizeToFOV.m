function [sensor,actualFOV] = sensorSetSizeToFOV(sensor,newFOV,scene,oi)
% Adjust sensor rows and columns so that horizontal FOV is deg (angle)
%
%   [sensor,actualFOV] = sensorSetSizeToFOV([sensor],newFOV,[scene],[oi])
%
% The aspect ratio of the sensor is left approximately unchanged. The
% imperfection is caused by rounding and the fact that the final size must
% be a multiple of the cfa size.
%
% The FOV of the sensor depends on the focal length to the optics and the
% size of the sensor. Hence, we normally send in the oi
% 
% We try to handle the human cone array case, which is special, by catching
% the string 'human' in the name and saying the block size is one.  This is
% managed in the sensorSet/Get operations.
%
%Example:
%  scene = sceneCreate; oi = oiCreate;
%  sensor = vcGetObject('sensor');
%  sensor = sensorSetSizeToFOV(sensor,1,scene,oi); % 1 deg field of view
%  [val,sensor] = vcGetSelectedObject('sensor'); 
%  sensor = sensorSetSizeToFOV(sensor,30,scene,oi); 
%  vcReplaceObject(sensor,val);
%
%  [sensor,actualFOV] = sensorSetSizeToFOV(sensor,3,scene,oi);
%
% Copyright ImagEval Consultants, LLC, 2005.

% PROGRAMMING:  Handle etendue.  Not addressed yet, sigh.

if notDefined('sensor'), sensor = vcGetObject('sensor'); end
if notDefined('newFOV'), error('horizontal field of view required'); end
if notDefined('scene'), scene = [];  end
if notDefined('oi'), oi = [];  end

% Get the size.  If it is 0, set to a small size.
sz = sensorGet(sensor,'size');
if isempty(sz)
    sz = [32,32];
    sensor = sensorSet(sensor,'size',sz);
end

% To compute the horizontal FOV, we need to know the distance from the
% sensor to the optics. Hence, we need to know oi and scene.
% If scene and oi are empty, then sensorGet uses the currently selected
% ones. If none are selected, then it uses some arbitrary default values.
% See the code in sensorGet.
currentFOV  = sensorGet(sensor,'fov horizontal',scene,oi);
newSize = round(sz * (newFOV/currentFOV) );

% The new sensor has to have at least the number of pixels in the cfa block
% pattern, and it has to be a multiple of the number of pixels in the block
% pattern. We make it slightly larger than absolutely necessary.

%%
cfaSize = sensorGet(sensor,'cfaSize');
if cfaSize ~= sz
    newSize = ceil(newSize ./ cfaSize) .* cfaSize;
    % If ceil(sz/cfaSize) is zero, we set size to one pixel cfa
    if newSize(1) == 0, newSize = cfaSize; end
end

% Set the new sizes
% For human case, the cone mosaic will be regenerated
sensor = sensorSet(sensor, 'size', newSize);
sensor = sensorClearData(sensor);

if nargout == 2, actualFOV = sensorGet(sensor,'fov'); end

end

