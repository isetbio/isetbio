function wvf = wvfComputePSF(wvf, showBar)
% Compute the psf for the wvf object. 
%
% Syntax:
%   wvf = wvfComputePSF(wvf, [showbar])
%
% Description:
%    If the psf is already computed and not stale, this will return fast.
%    Otherwise it computes and stores.
%
%    The point spread function is computed for each of the wavelengths
%    listed in the input wvf structure. The PSF computation is based on 10
%    orders of Zernike coefficients specified to the OSA standard.
%
%    The calculation also assumes that there is chromatic aberration of the
%    human eye, as embedded in the function wvfLCAFromWavelengthDifference, 
%    within the code in wvfComputePupilFunction.
%
%    Based on code provided by Heidi Hofer.
%
% Inputs:
%    wvf     - wavefront object
%    showbar - (Optional) Boolean dictating whether or not to show the
%              calculation wait bar
%
% Outputs:
%    wvf     - The wavefront object
%
% See Also:
%    wvfGet, wvfCreate, wvfSet, wvfComputePupilFunction, 
%    wvfLCAFromWavelengthDifference
%

% History:
%    08/20/11  dhb  Rename function and pull out of supplied routine.
%                   Reformat comments.
%    09/05/11  dhb  Rename. Rewrite for wvf i/o.
%    xx/xx/12       Copyright Wavefront Toolbox Team 2012
%    06/02/12  dhb  Simplify greatly given new get/set conventions.
%    07/01/12   bw  Adjusted for new wavelength convention
%    11/08/17  jnm  Comments & formatting
%

% Examples:
%{
    wvf = wvfCreate;
    wvf = wvfComputePSF(wvf)
%}

if notDefined('showBar'), showBar = false; end

% Only calculate if we need to -- PSF might already be computed and stored
if (~isfield(wvf, 'psf') || ~isfield(wvf, 'PSF_STALE') || ...
        wvf.PSF_STALE || ~isfield(wvf, 'pupilfunc') || ...
        ~isfield(wvf, 'PUPILFUNCTION_STALE') || wvf.PUPILFUNCTION_STALE) 
  
    % Initialize parameters. These are calc wave.
    wList = wvfGet(wvf, 'calc wave');
    nWave = wvfGet(wvf, 'calc nwave');
    pupilfunc = cell(nWave, 1);

    % Make sure pupil function is computed. This function incorporates the
    % chromatic aberration of the human eye.
    wvf = wvfComputePupilFunction(wvf, showBar);
    
    % wave = wvfGet(wvf, 'wave');
    psf = cell(nWave, 1);
    for wl = 1:nWave
        % Convert the pupil function to the PSF.
        % Requires only an ff2. 
        % Scale so that psf sums to unity.
        pupilfunc{wl} = wvfGet(wvf, 'pupil function', wList(wl));
        amp = fft2(pupilfunc{wl});
        inten = (amp .* conj(amp));   %intensity
        psf{wl} = real(fftshift(inten));
        psf{wl} = psf{wl} / sum(sum(psf{wl}));
    end
    
    wvf.psf = psf;
    wvf.PSF_STALE = false;
end

end
