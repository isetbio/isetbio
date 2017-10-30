function [temp, uv] = spd2cct(wave, spds, units)
% Convert a spectral power distribution to a correlated color temperature 
%
% Syntax:
%   [CCT, uv] = spd2cct(wave, spd, units)
%
% Description:
%    Calculates the correlated color temperature of a light from its
%    spectral power distribution in energy
%    CCT: Correlated color temperature.
%
% Inputs:
%	 wave  - Wavelengths of SPD.
%    spds  - Spectral power disbution of the lights. Can be in the columns
%            of a matrix.
%    units - SPD's units, to be converted to degrees Kelvin [Deprecated?]
%
% Outputs:
%    temp  - Color Temperature, in degrees Kelvin
%    uv    - Chromacity coordinates
%
% Notes:
%    * Are units still necessary? They aren't used anywhere in the function
%      or returned. 
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting

% Examples:
%{
   d = blackbody(400:10:700, 3500);
   spd2cct(400:10:700,d)

   d = blackbody(400:10:700, 6500);
   spd2cct(400:10:700,d)

   d = blackbody(400:10:700, 8500);
   spd2cct(400:10:700,d)
%}
XYZ = ieXYZFromEnergy(spds',wave);

% ISET returns uprime and vprime, which were defined in the 1960s. The flag
% makes sure we get 'uv' instead.
[u,v] =  xyz2uv(XYZ, 'uv');
uv = [u,v]';   % Format Jeff wrote for cct. u in first row, v in second
temp = cct( uv );

end
