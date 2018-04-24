function d = ieRad2deg(r, units)
% Convert radians to degrees
%
% Syntax:
%   d = ieRad2deg(r, units)
%
% Description:
%    Convert a measurement in Radians to Degrees. If 'arcmin' or 'arcsec'
%    are provided, convert to minutes.
%
% Inputs:
%    r     - Radians
%    units - (Optional) units, 'arcmin' or 'arcsec'
%
% Outputs:
%    d     - Calculated degrees
%
% Notes:
%    * [Note - DHB: Could consider deprecating as per the discussion in
%      issue 277.]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting
%    12/22/17  BW   Created, fixed error with notDefined, removed
%                   rad2deg
% Examples:
%{
    r = pi;
    ieRad2deg(r)
    r = ieDeg2rad(1);
    ieRad2deg(r, 'arcmin')
    ieRad2deg(r, 'arcsec')
%}

if notDefined('units')
    d = (180 / pi) * r; 
else
    d = r * ieUnitScaleFactor(units);
end   

end