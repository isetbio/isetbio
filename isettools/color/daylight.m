function [spec, XYZ] = daylight(wave, cct, units)
% Generate a daylight SPD with a correlated color temperature
%
% Syntax:
%   [spd, xyz] = daylight(WAVE, cct, units)
%
% Description:
%    Generates a daylight/sun spectral power distribution based on a
%    correlated color temperature (cct). 
%
% Inputs:
%    wave  - The wavelength vector in nm
%    cct   - The correlated color temperature
%    units - The corresponding units for the equation. Options are
%            'photons' or 'energy'
%
% Outputs:
%    spec  - The daylight spectral power distribution
%    xyz   - CIE XYZ values
%
% See Also:
%    blackbody, cct2sun
%

% History:
%    xx/xx/10       Copyright Imageval
%    10/30/17  jnm  Comments & formatting

% Examples
%{
   w = 400:700; 
   spd = daylight(w,6500,'energy');
   plot(w,spd)
%}
%{
   w = 400:700;
   spd = daylight(w,6500,'photons');
   plot(w,spd)
%}
%{
   w = 400:700;
   [spd, XYZ] = daylight(w,[4000 6500],'photons');
   plot(w,spd)
%}
%{
   wave = 400:770;         % Wavelength in nanometers
   cct = 4000:1000:10000;  % Correlated color temperature
   spd = daylight(wave,cct,'photons');
%}
if notDefined('wave')
    wave = 400:10:700;
end
if notDefined('units')
    units = 'energy';
end
if notDefined('cct')
    cct = 6500;
end

spec = cct2sun(wave, cct, units);

% Scale so first spectrum is 100 cd/m^2.
units = ieParamFormat(units);
switch units
    case 'photons'
        L = ieLuminanceFromPhotons(spec(:, 1)', wave);        
    case 'energy'
        L = ieLuminanceFromEnergy(spec(:, 1)', wave);
end
spec = (spec/L)*100;

if nargout == 2
    switch units
        case {'photons', 'quanta'}
            XYZ = ieXYZFromPhotons(spec', wave);
        case 'energy'
            XYZ = ieXYZFromEnergy(spec', wave);
    end
end

end