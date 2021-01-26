function [spotPattern] = drawSpot(spatialParams)
% Make an image with indicator values for a spot.
%
% Synopsis:
%    [spotPattern] = drawSpot(spatialParams)
%
% Description:
%    This function will draw a grayscale spot (increment) on top of a
%    background. Useful for simulating the spatial summation experiments of Davila and
%    Geisler (1991).
%
% Far part of image is set to 0, background is set to 1, spot to 2 in the returned image.
%
% The input should be a spatialParams structure.  See spatialParamsGenerate in the ISETBioCSF
% respository.

% History:
%   09/23/16    wst Wrote it
%   12/13/20    dhb Moved into isetbio.

% Set image values
backgroundValue = 1;
spotValue = 2;

% Check parameters
if (~strcmp(spatialParams.type,'Spatial') | ~strcmp(spatialParams.spatialType,'spot'))
    error('Wrong type of parameter structure passed');
end
if (spatialParams.row ~= spatialParams.col)
    error('This routine only works for square images');
end

% Center pixel
centerRow = round(spatialParams.row/2);
centerCol = round(spatialParams.col/2);
% Get background parameters
backgroundSizePixels = round((spatialParams.row/spatialParams.fieldOfViewDegs)*spatialParams.backgroundSizeDegs);

% Get spot parameters
spotRadiusDeg = spatialParams.spotSizeDegs/2;
spotRadiusPixels = (spatialParams.row/spatialParams.fieldOfViewDegs).*spotRadiusDeg;

% Create canvas in which to place the spot
spotPattern = zeros(spatialParams.row, spatialParams.col);

% Create and put background into canvas
backgroundImage = backgroundValue*ones(backgroundSizePixels,backgroundSizePixels);
if isodd(size(backgroundImage,1)) == 0
    minusHW = size(backgroundImage,1)/2-1;
    plusHW = (size(backgroundImage,1)/2);
else 
    minusHW = (size(backgroundImage,1)-1)/2;
    plusHW = minusHW;
end
spotPattern(centerRow-minusHW:centerRow+plusHW, centerCol-minusHW:centerCol+plusHW) = backgroundImage;

% Create spot in canvas.  Circle returns a tightly cropped boolean,
% with true in the spot pixels
croppedSpot = Circle(spotRadiusPixels);

% Pop spot into the image, preserving previous values
%
% Handle even versus odd pixel sizes
if isodd(size(croppedSpot,1)) == 0
    minusHW = size(croppedSpot,1)/2;
    plusHW = (size(croppedSpot,1)/2)-1;
else 
    minusHW = (size(croppedSpot,1)-1)/2;
    plusHW = minusHW;
end

temp = spotPattern(centerRow-minusHW:centerRow+plusHW, centerCol-minusHW:centerCol+plusHW);
temp(croppedSpot) = spotValue;
spotPattern(centerRow-minusHW:centerRow+plusHW, centerCol-minusHW:centerCol+plusHW) = temp;

end