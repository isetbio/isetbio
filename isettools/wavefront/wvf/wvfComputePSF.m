function wvf = wvfComputePSF(wvf, showBar)
% Compute the psf for the wvf object. 
%
%   wvf = wvfComputePSF(wvf)
%
% If the psf is already computed and not stale, this will return fast.
% Otherwise it computes and stores.
%
% The point spread function is computed for each of the wavelengths listed
% in the input wvf structure. The PSF computation is based on 10 orders of
% Zernike coefficients specified to the OSA standard.
%
% See also wvfGet, wvfCreate, wvfSet, wvfComputePupilFunction
%
% Based on code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.
% 9/5/11  dhb      Rename.  Rewrite for wvf i/o.
% 6/2/12  dhb      Simplify greatly given new get/set conventions.
% 7/1/12  bw       Adjusted for new wavelength convention
%
% Copyright Wavefront Toolbox Team 2012

if notDefined('showBar'), showBar = true; end

% Only calculate if we need to -- PSF might already be computed and stored
if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || ...
        wvf.PSF_STALE || ~isfield(wvf,'pupilfunc') || ...
        ~isfield(wvf,'PUPILFUNCTION_STALE') || wvf.PUPILFUNCTION_STALE) 
  
    % Initialize parameters.  These are calc wave.
    wList = wvfGet(wvf,'calc wave');
    nWave = wvfGet(wvf,'calc nwave');
    pupilfunc = cell(nWave,1);

    % Make sure pupil function is computed.  
    wvf = wvfComputePupilFunction(wvf, showBar);
    
    % wave = wvfGet(wvf,'wave');
    psf = cell(nWave,1);
    for wl = 1:nWave
        % Convert the pupil function to the PSF.
        % Requires only an ff2. 
        % Scale so that psf sums to unity.
        pupilfunc{wl} = wvfGet(wvf,'pupil function',wList(wl));
        amp = fft2(pupilfunc{wl});
        inten = (amp .* conj(amp));   %intensity
        psf{wl} = real(fftshift(inten));
        psf{wl} = psf{wl}/sum(sum(psf{wl}));
    end
    
    wvf.psf = psf;
    wvf.PSF_STALE = false;
end

end


