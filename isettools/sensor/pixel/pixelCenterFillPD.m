function isa = pixelCenterFillPD(isa,fillfactor)
% Adjust the pixel photodiode to be centered (fill factor [0,1])
%
%    isa = pixelCenterFillPD(isa,fillfactor)
%
% Create a centered photodetector with a particular fill factor within a
% pixel. (The pixel is attached to the sensor.)
%
% Example:
%  isa = pixelCenterFillPD(isa,0.5)
%
% See also: pixelPositionPD
%
% Copyright ImagEval Consultants, LLC, 2003.

% Programming notes:
%   This function should probably be written for the pixel, not the isa.
%
%      pixel = pixelCenterFillPD(pixel,fillfactor);
%
%  Also, changes in the size of the pixel within the GUI preserve the
%  fillfactor.  

if notDefined('isa'), [~,isa] = vcGetSelectedObject('ISA'); end
if notDefined('fillfactor'), fillfactor = 1;
elseif (fillfactor > 1) || (fillfactor < 0), 
    error('Fill factor must be between 0 and 1.  Parameter value = %f\n',fillfactor); 
end

pixel = sensorGet(isa,'pixel');

% Adjust pixel photodetector position to center with specified fill factor
% We define the fill factor as being the proportion of photodetector within
% the pixel plus the gap, not just the pixel.  
pixel.pdWidth  = sqrt(fillfactor)*pixelGet(pixel,'deltax');
pixel.pdHeight = sqrt(fillfactor)*pixelGet(pixel,'deltay');
isa = sensorSet(isa,'pixel',pixelPositionPD(pixel,'center'));

end