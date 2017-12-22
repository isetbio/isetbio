function d = rad2deg(r, units)
% Deprecated.  Use ieRad2deg or matlab version rad2deg.
%
% Syntax:
%   d = rad2deg(r, units)
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
%    * [Note: JNM - Changed varargin to units to support the function.
%      Fixed examples]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting

% Examples:
%{
    r = pi;
    rad2deg(r);
    rad2deg(r, "arcmin");
    rad2deg(r, "arcsec");
%}

error('Use ieRad2deg or matlab version of rad2deg');
end

% if isempty(units)
%     d = (180 / pi) * r; 
% else
%     d = r * ieUnitScaleFactor(units{1});
% end   
% 
% end