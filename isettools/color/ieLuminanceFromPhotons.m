function lum = ieLuminanceFromPhotons(photons,wave)
% Converts photons into energy and then calls ieLuminanceFromEnergy.
%
% Syntax:
%   lum = ieLuminanceFromPhotons(photons, wave)
%
% Description:
%    Calculate luminance (cd/m2) and related quantities (lux, lumens, cd)
%    from spectral power distribution in photon units. Which you get
%    depends on what units the spectral power distribution is in.
%        Luminance:          cd/m2 from photons/[sec-sr-m2-nm]
%        Illuminance:        lux from  photons/[sec-m2-nm]
%        Luminous flux:      lumens from photons/[sec-nm]
%        Luminous intensity: cd from photons/[sec-sr-nm]
%
%    This function contains examples of usage inline. To access these, type
%    'edit ieLuminanceFromPhotons.m' into the Command Window.
%
% Inputs:
%    photons - Matrix. The spectral power Distribution you wish to find the
%              luminance of.
%    wave    - Vector. The wavelengths, in nanometers
%
% Outputs:
%    lum     - Vector. The luminance, in candelas per meter squared, cd/m2
%              if input was in photons/[sec-sr-mr-nm]. See description.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   ieLuminanceFromEnergy
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/11/19  JNM  Formatting update

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