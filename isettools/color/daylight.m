function [spd, xyz] = daylight(wave, cct, units)
% Generate a daylight SPD with a correlated color temperature
%
% Syntax:
%   [spd, xyz] = daylight([wave], [cct], [units])
%
% Description:
%    Generates a daylight/sun spectral power distribution based on a
%    correlated color temperature (cct). 
%
%    All the returned spectra are normalized so that the luminance of the
%    first one is 100 cd/m2.
%
% Inputs:
%    wave  - The wavelength vector in nm (default 400:10:700).
%    cct   - The correlated color temperature (default 6500).
%    units - The corresponding units for the equation. Options are
%            'photons' or 'energy' (default 'energy').
%
% Outputs:
%    spd   - The daylight spectral power distribution, in columns of a
%            returned matrix.
%    xyz   - CIE XYZ values, in rows of returned matrix.
%
% Notes:
%   * [NOTE - DHB: Since spectra have been normalized so that first has a 
%      luminance of 100 cd/m2, should probably figure out what the
%      corresponding units of the spectra are and tell the user under
%      description above.
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
   [spd, xyz] = daylight(w,[4000 6500],'photons');
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
spd = (spd/L)*100;

if nargout == 2
    switch units
        case {'photons', 'quanta'}
            xyz = ieXYZFromPhotons(spd', wave);
        case 'energy'
            xyz = ieXYZFromEnergy(spd', wave);
        otherwise
            error('Unknown units specified');
    end
end

end