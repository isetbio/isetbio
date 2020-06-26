function uv = xyy2uv(xyy, format)
% Convert data from CIE xyY to CIE uv values
%
% Syntax:
%   uv = xyy2uv(xyy, [format])
%
% Description:
%    This routine calls xyy2xyz and then xyz2uv to convert from the
%    xyY representation to the uv representation.
%
%    Note that there are two closely related uv representations. The
%    current standard is called uprime, vprime and this is an
%    improvement made in the 1950s. This is what is returned by
%    xyz2uv by default.
%
%    This function contains examples of usage inline. To access these, type
%    'edit xyy2uv.m' into the Command Window.
%
% Inputs:
%    xyy    - Matrix. Nx3 [x, y, Y], chromacity Coordinates and Luminance.
%    format - String. Specifies the uv space. Default ''. Specifying the
%             default returns u-prime v-prime, but format = 'uv' returns
%             straight uv values. N.B. The 'uv' format is used for
%             correlated color temperature.
%
% Outputs:
%    uv     - Matrix. Nx2 [u-prime, v-prime] representation, or, if format
%             is 'uv', returns [u, v].
%
% Optional key/value pairs:
%    None.
%
% Examples are provided in the source code.
%
% See Also:
%   xyy2xyz, xyz2uv,cct2sun
%

% History
%	 12/12/17  BW   first draft
%    07/15/19  JNM  Formatting

% Examples:
%{
    xyY = [.3221 .3322 100];
    xyy2uv(xyY, 'uv')  % u, v used for cct2sun
    xyy2uv(xyY)        % uprime, vprime format
%}

% Default is uprime, vprime format
if notDefined('format'), format = ''; end

% Check input
if size(xyy, 2) ~= 3
    error('Input must be Nx3, with [x, y, Y] in the rows.')
end

XYZ = xyy2xyz(xyy);
uv = xyz2uv(XYZ, format);  % Set 'uv' flag for cct2sun or cct case.

end
