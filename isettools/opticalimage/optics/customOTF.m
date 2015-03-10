function [OTF2D, fSupport] = customOTF(oi,fSupport,wavelength,units)
%Interpolate optics OTF for shift-invariant calculation in optical image
% 
%  [OTF2D,fSupport] = customOTF(oi,[fSupport],[wavelength = :],[units='mm'])
%
% In the shift-invariant optics model, custom data are stored in the
% optics.OTF slot.  This routine reads the otf data and interpolates them
% to the fSupport and wavelength of the optical image or optics structure.
%
% The returned OTF is normalized so that all energy is transmitted (i.e.,
% the DC value is 1).  This is done by normalizing the peak value to one.
% If we ever have a case when the peak is other than the DC, we have a
% problem with energy conservation - where did the photons go?
%
% The units for the frequency support are cycles/millimeters.  
% Perhaps we should add a 'units' input argument here. 
% 
% Examples:
%
% See also: oiCalculateOTF, dlCore.m
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('oi'),         error('Optical image required.'); end
if notDefined('wavelength'), wavelength = oiGet(oi,'wavelength'); end

% In the custom case, I think the units should always be millimeters.
if notDefined('units'),      units = 'mm'; end  
if notDefined('fSupport'),   fSupport = oiGet(oi,'fSupport',units); end

fx    = fSupport(:,:,1); fy = fSupport(:,:,2);
nX    = size(fx,2);      nY = size(fy,1);
nWave = length(wavelength);

optics     = oiGet(oi,'optics');

% This is the representation of the OTF
% The units should be specified and optional; we are using mm for
% now.  We can run into trouble if the values of fx or fy exceed those in X
% and Y.  This produces NaNs in the interpolated OTF2D, below.  We could
% check conditions on these variables:
%   max(abs(X(:))), max(abs(fx(:)))
%   max(abs(fy(:))), max(abs(Y(:)))
% I would like to replace the opticsGet here with an oiGet call with
% frequencySupport.  There are some inconsistencies that I don't understand
% and we have to clarify.
% BW, 2010 May
otfSupport = opticsGet(optics,'otfSupport');  
[X, Y]     = meshgrid(otfSupport{1},otfSupport{2});

% set interpolation method
if ieSessionGet('gpu compute')
    method = 'linear';
    fx = gpuArray(fx); fy = gpuArray(fy);
else
    method = '*linear';
end

% Find the OTF at each wavelength. We may be interpolating from the custom
% data.
if isscalar(wavelength)
    % Should we be interpolating here?
    OTF2D = opticsGet(optics, 'otfData', wavelength);
    % figure(1); mesh(X,Y,OTF2D); 
    % figure(1); mesh(X,Y,abs(OTF2D))
    
    % Do interpolation, use fftshift to take care of DC positions
    OTF2D = ifftshift(interp2(X, Y, fftshift(OTF2D), fx, fy, method,0));

else
    % Do it wavelength by wavelength
    if ieSessionGet('gpu compute')
        OTF2D = zeros(nY, nX, nWave, 'gpuArray');
    else
        OTF2D = zeros(nY,nX,nWave);
    end
    for ii=1:length(wavelength)
        tmp = opticsGet(optics, 'otfData', wavelength(ii));
        OTF2D(:,:,ii) = ...
            fftshift(interp2(X, Y, fftshift(tmp), fx, fy, method,0));
    end
end

end