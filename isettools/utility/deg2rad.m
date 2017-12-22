function r = deg2rad(d, units)
% renamed to ieDeg2rad and depreciated because of conflict with Matlab
% deg2rad. 
%
% Convert degrees to radians
%
% Syntax:
%   r = deg2rad(d, [units])
%
% Description:
%    Convert degrees to radians. Apply additional conversions if
%    appropriate units are provided. (Ex. arcmin, arcsec)
%
% Inputs:
%    d     - Degrees to convert
%    units - Units to apply to radians in conversion. 
%
% Outputs:
%    r     - measurement in Radians (unless other unit specified)
%
% Notes:
%    * [Note: XXX - Perhaps this routine should also adopt the varargin
%      structure in that function.]
%    * [Note: JNM - Added in units section (that changed from rad2deg) -
%      please check over my work!]
%
% See Also:
%    rad2deg

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/01/17  jnm  Formatting, Add units argument, if statement based on
%                   units argument, and scale multiplier. 

% Examples:
%{
    deg2rad(90)
    deg2rad(90, 'arcmin')
%}
error('Use ieDeg2rad or matlab version of deg2rad');
end

% if notDefined('units')
%     r = (pi / 180) * d;
% else
%     r = (pi / 180) * d * ieUnitScaleFactor(units);
% end
