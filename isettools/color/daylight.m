function [spd, XYZ] = daylight(wave, cct, units)
% Generate a daylight SPD with a correlated color temperature (cct)
%
% Syntax:
%   [spd, XYZ] = daylight([wave], [cct], [units])
%
% Description:
%    Generates a daylight spectral power distribution based on a
%    correlated color temperature (cct).
%
%    Tthe returned spectra are normalized so that the luminance of the
%    first returned spd is 100 cd/m2.  See the 3rd example below for how to
%    calculate the luminance of each of the lights.
%
% Inputs:
%    wave  - The wavelength vector in nm (Default 400:10:700).
%    cct   - The correlated color temperature (Default 6500).
%    units - The corresponding units for the equation. Options are
%            'photons' or 'energy' (Default 'energy').
%
% Outputs:
%    spd   - The daylight spectral power distribution, in columns of a
%            returned matrix.
%    XYZ   - CIE XYZ values, in rows of returned matrix.
%
% See Also:
%    blackbody, cct2sun
%

% History:
%    xx/xx/10       Copyright Imageval
%    10/30/17  jnm  Comments & formatting
%    11/13/17  baw  Added example for computing luminance of all spds
%    11/16/17  jnm  Formatting
%

% Examples
%{
   w = 400:700; 
   spd = daylight(w, 6500, 'energy');
   plot(w, spd)
%}
%{
   w = 400:700;
   spd = daylight(w, 6500, 'photons');
   plot(w, spd)
%}
%{
   w = 400:700;
   [spd, xyz] = daylight(w, [4000:1000:8000], 'photons');
   plot(w, spd)
   % luminance of all the lights
   L = ieLuminanceFromPhotons(spd', w);   
%}
%{
   wave = 400:770;         % Wavelength in nanometers
   cct = 4000:1000:10000;  % Correlated color temperature
   spd = daylight(wave, cct, 'photons');
   plot(w, spd)
%}

if notDefined('wave'), wave = 400:10:700; end
if notDefined('units'), units = 'energy'; end
if notDefined('cct'), cct = 6500; end

% Make the spd
spd = cct2sun(wave, cct, units);

% Scale so first spectrum is 100 cd/m^2.
units = ieParamFormat(units);
switch units
    case 'photons'
        L = ieLuminanceFromPhotons(spd(:, 1)', wave);        
    case 'energy'
        L = ieLuminanceFromEnergy(spd(:, 1)', wave);
end
spd = (spd / L) * 100;

if nargout == 2
    switch units
        case {'photons', 'quanta'}
            XYZ = ieXYZFromPhotons(spd', wave);
        case 'energy'
            XYZ = ieXYZFromEnergy(spd', wave);
        otherwise
            error('Unknown units specified');
    end
end

end