function oi = wvf2oi(wvf)
% Convert wavefront data to ISETBIO optical image with optics
%
% Syntax:
%   oi = wvf2oi(wvf)
%
% Description:
%    Use Zernicke polynomial data in the wvfP structure and create an
%    ISETBIO optical image whose optics match the wavefront data structure.
%
%    Before calling this function, compute the PSF of the wvf structure.
%
%    Non-optics aspects of the oi structure take on default values.
%
% Inputs:
%    wvf - A wavefront parameters structure (with a computed PSF)
%
% Outputs:
%    oi  - ISETBIO optical image
%
% See Also:
%    oiCreate, oiPlot
%

% History:
%	 xx/xx/12       Copyright Wavefront Toolbox Team 2012
%    11/13/17  jnm  Comments & formatting
%    01/01/18  dhb  Set name and oi wavelength from wvf.
%              dhb  Check for need to interpolate, skip if not.

% Examples
%{
    wvf = wvfCreate;
    wvf = wvfComputePSF(wvf);
    oi = wvf2oi(wvf);
    oiPlot(oi, 'psf550');
%}
%{
    wvf = wvfCreate('wave',[400 550 700]');
    wvf = wvfSet(wvf, 'zcoeff', 1, 'defocus');
	wvf = wvfComputePSF(wvf);
    oi = wvf2oi(wvf);
    oiPlot(oi, 'psf550');
%}

%% Set up parameters
if notDefined('wvf'), error('Wavefront structure required.'); end
wave = wvfGet(wvf, 'calc wave');

%% First we figure out the frequency support.
fMax = 0;
for ww=1:length(wave)
    f = wvfGet(wvf, 'otf support', 'mm', wave(ww));
    if max(f(:)) > fMax
       fMax = max(f(:));
       maxWave = wave(ww);
    end
end

% Make the frequency support in ISET as the same number of samples with the
% wavelength with the highest frequency support from WVF.
%
% This support is set up with sf 0 at the center of the returned vector,
% which matches how the wvf object returns the otf.
%
% This section is here in case the frequency support for the wvf otf varies
% with wavelength.  Not sure that it ever does. There is a conditional that
% skips the inerpolation if the frequency support at a wavelength matches
% that with the maximum, so this doesn't cost us much time.
fx = wvfGet(wvf, 'otf support', 'mm', maxWave);
fy = fx;
[X, Y] = meshgrid(fx, fy);
c0 = find(X(1, :) == 0);
tmpN = length(fx);
if (floor(tmpN/2)+1 ~= c0)
    error('We do not understand where sf 0 should be in the sf array');
end

%% Set up the OTF variable for use in the ISETBIO representation
nWave = length(wave);
nSamps = length(fx);
otf = zeros(nSamps, nSamps, nWave);

%% Interpolate the WVF OTF data into the ISET OTF data for each wavelength.
%
% The interpolation seems to be here in case there is different frequency
% support in the wvf structure at different wavelengths.
for ww=1:length(wave)
    f = wvfGet(wvf, 'otf support', 'mm', wave(ww));
    if (f(floor(length(f)/2)+1) ~= 0)
        error('wvf otf support does not have 0 sf in the expected location');
    end
    
    % Apply fftshift to convert otf to DC in center, so interp will work right.
    thisOTF = fftshift(wvfGet(wvf, 'otf', wave(ww)));
    if (all(f == fx))
        est = thisOTF;
    else
        est = interp2(f, f', thisOTF, X, Y, 'cubic', 0);
    end
    
    % Isetbio wants the otf with (0,0) sf at the upper left.  We
    % accomplish this by applying ifftshift to the wvf centered format.
    wvf2oiBackCompat = false;
    if (ispref('isetbioBackCompat','wvf2oi'))
        if (getpref('isetbioBackCompat','wvf2oi'))
            wvf2oiBackCompat = true;
        end
    end
    if (wvf2oiBackCompat)
        % This is the old way.
        %
        % It is tragic that fftshift does not shift so that the DC term is in
        % (1, 1). Rather, fftshift puts the DC at the the highest position.
        % So, we don't use this
        %
        %   otf(:, :, ww) = fftshift(otf(:, :, ww));
        %
        % Rather, we use circshift. This is also the process followed in the
        % psf2otf and otf2psf functions in the image processing toolbox. Makes
        % me think that Mathworks had the same issue. Very annoying. (BW)
        
        % We identified the (r, c) that represent frequencies of 0 (i.e., DC).
        % We circularly shift so that that (r, c) is at the (1, 1) position.
        r0 = c0;
        otf(:, :, ww) = circshift(est, -1 * [r0 - 1, c0 - 1]);
    else
        % This is the new way, which seems cleaner to me (DHB).
        otf(:, :, ww) = ifftshift(est);
    end
end

%% Place the frequency support and OTF data into an ISET structure.
%
% Build template with standard defaults
oi = oiCreate;
oi = oiSet(oi,'name',wvfGet(wvf,'name'));

% Copy the OTF parameters.
oi = oiSet(oi, 'optics OTF fx', fx);
oi = oiSet(oi, 'optics OTF fy', fy);
oi = oiSet(oi, 'optics otfdata', otf);
oi = oiSet(oi, 'optics OTF wave', wave);
oi = oiSet(oi, 'wave', wave);

end
