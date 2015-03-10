function oi = oiPad(oi,padSize,sDist,direction)
%Pad the oi irradiance data with zeros, usually prior to applying OTF
%
%     oi = oiPad(oi,padSize,[sDist],direction)
%
% For optics calculations we need to pad the size (to avoid edge wrapping).
% Here we pad the spatial dimensions with 0s. By changing the row and
% column numbers, we also must and adjust some parameters, such as the
% horizontal field of view accordingly.
%
% You can set the argument direction = 'both','pre', or 'post' to pad both
% or only on one side. By default, the zero-padding takes place on all
% sides of the image.  Thus,  by default if padSize(1) is 3, we add 3 rows
% on top and bottom (total of 6 rows).
%
% Example:
%   oi = oiPad(oi,[8,8,0]);
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('sDist'), 
    scene = vcGetObject('scene');
    if isempty(scene)
        warning('oiPad: No scene, assuming 1 m sDist');
        sDist = 1;
    else  sDist = sceneGet(scene,'distance'); 
    end
end
if notDefined('direction'), direction = 'both'; end

% We make sure padSize matches the dimensionality of photons.
% Probably not necessary.  But ...
if ismatrix(padSize), padSize(3) = 0; end


photons = oiGet(oi,'photons');
% To prevent ieCompressData error, we set the surrounding region as the
% mean of the data at each wavelength samples. In this way, we will have
% the mean luminance matched.
if isa(photons, 'gpuArray')
    padval = gather(mean(mean(photons)));
else
    padval = mean(mean(photons));
end


% Pad one wavelength at a time at single
% Compute the size of the new, padded array.
photons = single(photons);
[r,c] = size(padarray(photons(:,:,1), padSize, padval(1), direction));

% Figure out the number of wavebands
nwave = size(photons,3);

% Now, use single instead of double precision.
newPhotons = zeros(r, c , nwave, 'single');
for ii = 1 : nwave
    newPhotons(:,:,ii) = ...
        padarray(photons(:,:,ii), padSize, padval(ii), direction);
end

% Clear unused stuff
clear photons;
photons = newPhotons;
clear newPhotons;

% The sample spacing of the optical image at the surface of the sensor must
% be adjusted for the padding.  We must make this adjustment before putting
% the new data into the oi because we need to preserve the number of
% columns until we are done with this calculation.

% The width per horizontal sample at the sensor surface is the ratio of the
% width to the number of columns.  The new number of columns is the sum of
% the current number and the horizontal pad size, which is in pad(2)
if strcmp(direction, 'both'), padCols = padSize(2)*2;
else padCols = padSize(2); end
newWidth = oiGet(oi,'width')* (1 + padCols/oiGet(oi, 'cols'));

% Find the distance from the image to the lens
imageDistance = oiGet(oi, 'optics imageDistance', sDist);

% Now we compute the new horizontal field of view using the formula that
% says the opposite over adjacent is the tangent of the angle.  We return
% the value in degrees
oi = oiSet(oi, 'h fov', 2*atand(newWidth/2/imageDistance));

% Now we adjust the columns by placing in the new photons
oi = oiSet(oi, 'photons', photons);

end