function [xySubMosaic, coneSubMosaic, support, spread, delta] = ...
    conePlotHelper(sensor, conesToPlot, coneSize)
% Return the input variables used by the conePlot function
%
% [xySubMosaic, coneSubMosaic, support, spread, delta] = ...
%     conePlotHelper(sensor, conesToPlot, coneSize)
%
% Inputs:
%   Sensor
%   Desired number of cones from center to plot [default: max sensor size]
%   Cone diameter as image pixels [default: 10] 
%
% Returns:
%   xySubMosaic
%   coneSubMosaic
%   support
%   spread
%   delta
%   
% Example: 
%   Generate an image of the central 30x30 cones, with each cone having a
%   20 pixel diameter
%
%   sensor = sensorCreate('human');
% 
%   conesToPlot = [];
%   coneSize = 20;
% 
%   [xy, coneType, support, spread, delta] = conePlotHelper(sensor, conesToPlot, coneSize);
% 
%   whiteBackground = false;
%   [~,~,~, coneMosaicImage] = conePlot(xy,coneType, support,spread,delta, whiteBackground);
%   imshow(coneMosaicImage); truesize;
%
% See also: sensorPlot
%
% NC, BW, maybe others, ISETBIO Team, Copyright 2015

% We should be able to multiply this image by the cone absorptions to
% produce a nice colored image.
if notDefined('sensor'), error('sensor structure required.'); end
if notDefined('conesToPlot'), conesToPlot = sensorGet(sensor,'size'); end
if isscalar(conesToPlot), conesToPlot(2) = conesToPlot(1); end
if notDefined('coneSize'), coneSize = 10; end

separation(1) = sensorGet(sensor, 'wspatial resolution', 'um');
separation(2) = sensorGet(sensor, 'hspatial resolution', 'um');

xy = sensorGet(sensor, 'cone xy');
xCoords = squeeze(xy(:,1));
yCoords = squeeze(xy(:,2));

% These are the cones we will plot.  They are chosen from the center of the
% image.
selectConeIndices = find( ...
    (abs(xCoords) < separation(1)*conesToPlot(2)/2) & ...
    (abs(yCoords) < separation(2)*conesToPlot(1)/2) ...
    );

coneType = sensorGet(sensor, 'coneType');

spread        = 0.25*coneSize; 
support       = round(spread*4*[1 1]); 
delta         = 1.5/coneSize;
coneSubMosaic = coneType(selectConeIndices);
xySubMosaic   = xy(selectConeIndices,:);

end
