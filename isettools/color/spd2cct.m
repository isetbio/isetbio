function [CCT, uv] = spd2cct(wave, spds)
% Convert a spectral power distribution to a correlated color temperature 
%
% Syntax:
%   [CCT, uv] = spd2cct(wave, spd)
%
% Description:
%    Calculates the correlated color temperature of a light from its
%    spectral power distribution in energy units, result is in degrees
%    Kelvin.
%    CCT: Correlated color temperature.
%
% Inputs:
%	 wave  - Wavelengths of SPD.
%    spds  - Spectral power disbution of the lights. Can be in the columns
%            of a matrix.
%
% Outputs:
%    CCT   - Correlated color temperature, in degrees Kelvin
%    uv    - Chromacity coordinates
%
% See also:
%    cct, xyz2uv

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting
%    10/31/17  dhb  Remove unused units arg, as per note query.  Deleted
%                   the note.

% Examples:
%{
   d = blackbody(400:10:700, 3500);
   spd2cct(400:10:700,d)

   d = blackbody(400:10:700, 6500);
   spd2cct(400:10:700,d)

   d = blackbody(400:10:700, 8500);
   spd2cct(400:10:700,d)
%}

% Convert to XYZ
XYZ = ieXYZFromEnergy(spds',wave);

% ISET returns uprime and vprime, which were defined in the 1960s. The flag
% makes sure we get 'uv' instead.
[u,v] =  xyz2uv(XYZ, 'uv');

% Format for cct routine. u in first row, v in second
uv = [u,v]';  

% Get temp
CCT = cct( uv );

end
