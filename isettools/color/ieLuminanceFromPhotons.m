function lum = ieLuminanceFromPhotons(photons,wave)
% Converts photons into energy and then calls ieLuminanceFromEnergy.
%
% Syntax:
%   lum = ieLuminanceFromPhotons(photons, wave)
%
% Description:
%    Calculate luminance (cd/m2) and related quantities (lux, lumens, cd)
%    from spectral photons
%
% Inputs:
%    photons - the Spectral Power Distribution you wish to find the
%              luminance of.
%    wave    - The wavelengths, in nanometers
%
% Outputs:
%    lum     - The luminance, in candelas per meter squared, cd/m2
%
% See Also:
%   ieLuminanceFromEnergy


% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting

% Examples:
%{
   wave = 400:770;
   cct = 4000:1000:10000;
   spd = daylight(wave, cct, 'photons');

   % Calculate the luminance of these SPDs
   lum = ieLuminanceFromPhotons(spd', wave(:));
%}
energy = Quanta2Energy(wave, photons);
lum = ieLuminanceFromEnergy(energy, wave);

end