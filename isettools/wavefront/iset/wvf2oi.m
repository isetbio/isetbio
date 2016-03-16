function oi = wvf2oi(wvf)
% Convert wavefront data to ISETBIO optical image with optics
%
%    oi = wvf2oi(wvf)
%
% Use Zernicke polynomial data in the wvfP structure and create an ISETBIO
% optical image whose optics match the wavefront data structure.
%
% Before calling this function, compute the PSF of the wvf structure.
%
% wvf:  A wavefront parameters structure (with a computed PSF)
%
% Examples
%   wvf = wvfCreate; wvf = wvfComputePSF(wvf);
%   oi = wvf2oi(wvf); oiPlot(oi,'psf550');
%
%   wvf = wvfSet(wvf,'zcoeff',1,'defocus'); wvf = wvfComputePSF(wvf);
%   oi = wvf2oi(wvf); oiPlot(oi,'psf550');
%
% See also: 
%   oi = oiCreate('wvf human'); 
%   oiPlot(oi,'psf550'); oiPlot(oi,'psf',[],420);
%
% Copyright Wavefront Toolbox Team 2012

%% Set up parameters
if notDefined('wvf'), error('Wavefront structure required.'); end

wave = wvfGet(wvf,'calc wave');

% First we figure out the frequency support.
fMax = 0;
for ww=1:length(wave)
    f = wvfGet(wvf,'otf support','mm',wave(ww));
    if max(f(:)) > fMax
       fMax = max(f(:)); maxWave = wave(ww);
    end
end

% Make the frequency support in ISET as the same number of samples with the
% wavelength with the highest frequency support from WVF.
fx = wvfGet(wvf,'otf support','mm',maxWave);
fy = fx;
[X,Y] = meshgrid(fx,fy);
c0 = find(X(1,:) == 0); r0 = find(Y(:,1) == 0);

%% Set up the OTF variable for use in the ISETBIO representation

nWave = length(wave); nSamps = length(fx);
otf = zeros(nSamps,nSamps,nWave);

% Interpolate the WVF OTF data into the ISET OTF data for each wavelength.
for ww=1:length(wave)
    f = wvfGet(wvf,'otf support','mm',wave(ww));
    thisOTF = wvfGet(wvf,'otf',wave(ww));
    est = interp2(f,f',thisOTF,X,Y,'cubic',0);
    
    % It is tragic that fftshift does not shift so that the DC term is in
    % (1,1). Rather, fftshift puts the DC at the the highest position.
    % So, we don't use this
    %
    %   otf(:,:,ww) = fftshift(otf(:,:,ww));
    %
    % Rather, we use circshift.  This is also the process followed in the
    % psf2otf and otf2psf functions in the image processing toolbox.  Makes
    % me think that Mathworks had the same issue.  Very annoying. (BW)
    
    % We identified the (r,c) that represent frequencies of 0 (i.e., DC).
    % We circularly shift so that that (r,c) is at the (1,1) position.
    otf(:,:,ww) = circshift(est,-1*[r0-1,c0-1]);
    
end

% I sure wish this was real all the time. Sometimes (often?) it is. 
%  psf = otf2psf(otf(:,:,ww));
%  if ~isreal(psf), disp('psf not real'); end
%  vcNewGraphWin; mesh(psf)

%% Place the frequency support and OTF data into an ISET structure.

% Build template
oi = oiCreate;

% Same wavelengths
% oi = oiSet(oi,'wave',wave);

% Copy the OTF parameters.
% Note the awful repetition of setting the wave.  Sigh. This should get
% better. 
oi = oiSet(oi,'optics OTF fx', fx);
oi = oiSet(oi,'optics OTF fy', fy);
oi = oiSet(oi,'optics otfdata', otf);
% oi = oiSet(oi,'optics wave',wave);
oi = oiSet(oi,'optics OTF wave',wave);

% Z = wvfGet(wvf,'zcoeffs');
% oi = oiSet(oi,'zernike',Z);

end
