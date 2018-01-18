function r = ieDeg2rad(d, units)
% Convert degrees to radians
%
% Syntax:
%   r = ieDeg2rad(d, [units])
%
% Description:
%    Convert degrees to radians. Apply additional conversions if
%    appropriate units (arcmin, arcsec) are provided.
%
% Inputs:
%    d     - Degrees to convert
%    units - Units to apply to radians in conversion. 
%
% Outputs:
%    r     - measurement in Radians (unless other unit specified)
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: DHB - Could consider deprecating as per our discussion in
%      issue 277.]
%
% See Also:
%    rad2deg, deg2rad
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/01/17  jnm  Formatting, Add units argument, if statement based on
%                   units argument, and scale multiplier. 
%    12/22/17  BW   Introduced, removed deg2rad to avoid matlab
%                   conflict, replaced throughout code.
%    01/11/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    ieDeg2rad(90)
    ieDeg2rad(1, 'arcmin')
    ieDeg2rad(1, 'arcsec')
%}

if notDefined('units')
    r = (pi / 180) * d;
else
    r = (pi / 180) * d * ieUnitScaleFactor(units);
end
