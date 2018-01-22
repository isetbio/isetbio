function uv = xyy2uv(xyy,format)
% Convert data from CIE xyY to CIE uv values
%
% Syntax:
%   uv = xyy2uv(xyy,[format])
%
% Description:
%    This routine calls xyy2xyz and then xyz2uv to convert from the
%    xyY representation to the uv representation. 
%
%    Note that there are two closely related uv representations.  The
%    current standard is called uprime, vprime and this is an
%    improvement made in the 1950s.  This is what is returned by
%    xyz2uv by default.  
%
% Inputs:
%    xyy - Nx3 [x,y,Y], chromacity Coordinates and Luminance 
%    format - Specifies the uv space. By default format ='' returns
%             u-prime v-prime, but format = 'uv' returns striaght uv values.
%             N.B. The 'uv' format is used for correlated color temperature.
%
% Outputs:
%    uv - Nx2 [u-prime,v-prime] representation, or [u,v] if format='uv'
%
% Optional key/value pairs:
%    None.
%
% Examples are provided in the source code.
%
% See also: xyy2xyz, xyz2uv,cct2sun

% History
%    12.12.2017 BW first draft

% Examples:
%{
   xyY = [.3221 .3322 100];
   xyy2uv(xyY,'uv')   % u,v used for cct2sun
   xyy2uv(xyY)        % uprime, vprime format
%}

% Default is uprime, vprime format
if notDefined('format'),format = ''; end

% Check input
if size(xyy, 2) ~= 3
    error('Input must be Nx3, with x,y,Y in the rows.')
end

XYZ = xyy2xyz(xyy);

uv = xyz2uv(XYZ,format);   % Set 'uv' flag for cct2sun or cct case.

end
