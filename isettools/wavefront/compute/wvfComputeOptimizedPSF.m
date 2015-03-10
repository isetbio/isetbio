function wvfParams = wvfComputeOptimizedPSF(wvfParams)
% wvfParams = wvfComputeOptimizedPSF(wvfParams)
%
% Optimize the PSF seen at a specified wavelength.  Optimization is performed on
% the defocus parameter, relative to a specified nominal focus wavelength.
% The full polychromatic PSF is returned at the list of specified
% wavelengths.
%
% This is implemented as a call into wvfComputeOptimzedConePSF.
%
% Required input fields for wvfParams struct - see comment in
% wvfComputePupilFunction for more details.
%
%   criterionFraction - The figure of merit is the radius of a centered and circularly averaged version of
%                       the psf that contains the specified criterionFraction of the total mass of the PSF seen by each cone.
%                       The smaller this radius, the better.
%   optimizeWl -        Wavelenght to optimize for.
%   wls -               Column vector of wavelengths over which polychromatic psf is computed.
%   zcoeffs -           Zernike coefficients.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%   caclpupilsize -     Size over which returned pupil function is calculated, in MM.
%   wls -               Column vector of wavelengths over which to compute, in NANOMETERS.
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   fieldSampleSizeMMperPixel - Size in mm of each pixel of the pupil
%                       function.
%   sizeOfFieldMM -     Size of square image over which the pupile function is computed in MM.
%
% Optional input fields for wvfParams struct
%   sceParams -         Parameter structure for Stiles-Crawford correction.  If missing or set to empty,
%                       no correction and is set to empty on return.
%
% Output fields set in wvfParams struct
%   conepsf -           Calcuated psf for each cone in T_cones, third dimension indexes cone type.
%   defocusDiopters -   The defocus added in to optimize.
%   coneSceFrac -       Vector with calculated SCE fraction for each cone type.
%   psf -               Calcuated polychromatic psf. Third dimension of returned matrix indexes wavelength.
%   pupilfunc -         Calculated pupil function.  Third dimension of returned matrix indexes wavelength
%   arcminperpix -      Arc minutes per pixel for returned psfs.
%   strehl -            Strehl ratio of psf at each wavelength.  If SCE correction is specified, the returned
%                       strehl ratio is to the diffraction limited psf with the same SCE assumed.
%   sceFrac -           Fraction of light actually absorbed when SCE is taken into account, at each wavelength.
%   areapix -           Number of pixels within the computed pupil aperture at each wavelength
%   areapixapod -       Number of pixels within the computed pupil aperture at each wavelength,
%                       multiplied by the Stiles-Crawford aopdization.
%   defocusMicrons -    Defocus added in to zcoeffs(4) at each wavelength, in microns.
%
% 9/9/11   dhb  Wrote it.

index = find(wvfParams.optimizeWl == wvfParams.wls);
if (isempty(index))
    error('Desired wavelength to optimize for not inluded in wavelengths to compute for.');
end

% Set up the fields we need to make wvfComputeOptimizedConePSF do what we want.
wvfParams.coneWeights = 1;
wvfParams.T_cones = zeros(1,length(wvfParams.wls));
wvfParams.T_cones(index) = 1;
wvfParams.weightingSpectrum = zeros(length(wvfParams.wls),1);
wvfParams.weightingSpectrum(index) = 1;

% Do the work
wvfParams = wvfComputeOptimizedConePSF(wvfParams);

% Remove the fields we added in and no longer need
wvfParams = rmfield(wvfParams,'conepsf');
wvfParams = rmfield(wvfParams,'coneSceFrac');
wvfParams = rmfield(wvfParams,'coneWeights');
wvfParams = rmfield(wvfParams,'T_cones');
wvfParams = rmfield(wvfParams,'weightingSpectrum');