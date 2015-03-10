function [img,vci] = displayRender(img,vci,sensor)
% Transform rgb image data from the internal color space to display space
%
%   [img,vci] = displayRender(img,vci,sensor);
%
% The rendered display image RGB is always between 0 and 1. We use the
% linear display primary representation as the default RGB values.
%
% In the general processing pipeline, sensor data are demosaicked,
% converted to the internal color space, and then color balanced.  After
% the color balancing the data must be converted from the internal color
% space to the display representation.  It is that last step that is
% performed here.
% 
% The transform from the internal color space to the display values is
% computed using ieInternal2Display.  This transform is stored in the vci
% list of transforms. 
%
% The ratio of the image data maximum to the display absolute maximum is
% set equal to the ratio of the sensor data maximum to the sensor's
% absolute maximum. So, suppose that the the largest sensor value is 0.8
% volts and the voltage swing of the sensor is 1.5 volt. Then, if the
% maximum display output is 1 (as it always is), we set the maximum in this
% particular image to be (0.8/1.5). 
%
% The purpose of the scaling is to make sure that when we saturate the
% sensor, the output image maps into the peak of the display response.
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('img'), error('Must define image'); end
if notDefined('vci'), vci = vcGetObject('vcimage'); end
if notDefined('sensor'), sensor = vcGetObject('sensor'); end

% There are a couple of other ways we can do the rendering

% This stage converts from the internal color space to the display color
% space (linearly).
switch ieParamFormat(imageGet(vci,'internalCS'))
    case {'xyz','stockman'}
        % This is the transform for several calibrated color space
        M = ieInternal2Display(vci);
    case {'sensor'}
        % nSensors = imageGet(vci,'nSensorInputs');
        % No color space conversion in this case.
        % Simply copy the data from the sensor space to display RGB
        N = size(imageGet(vci,'illuminant correction transform'),2);
        M = eye(N,3);  % Always three display outputs (RGB).
    otherwise
        error('Unknown internal color space')
end

method = ieParamFormat(imageGet(vci,'Sensor conversion method'));
switch lower(method)
    case {'current','currentmatrix','manualmatrixentry','none'}
        % Do not apply another matrix if the user set the Transform manually
    case {'sensor','mccoptimized', 'esseroptimized'}
        % Apply the final transform from ICS to Display
        % M   = scaleMTransform(M,vci,sensor);
        vci = imageSet(vci,'ics2display',M);
        img = imageLinearTransform(img,M);
        
    otherwise
        error('Unknown color conversion method %s',method)
end

% The display image RGB is always between 0 and 1. Set the maximum
% image value to the ratio of the maximum value in this sensor
% image divided by the sensor's maximum possible value.
imgMax = max(img(:));
img = (img/imgMax)*sensorGet(sensor,'response ratio');

% figure; imagescRGB(vci.data.result);

end

        