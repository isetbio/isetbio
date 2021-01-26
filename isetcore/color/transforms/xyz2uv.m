function [u, v] = xyz2uv(xyz, format)
% Convert CIE XYZ to uv chromaticity coordinates
%
% Syntax:
%   [u, v] = xyz2uv(xyz, [format])
%
% Description:
%    Convert XYZ to uprime, vprime chromaticity coordinates (uniform
%    chromaticity space).
%
%    XYZ contains the values in the rows.
%    format: If you want the (u, v), and not (uprime, vprime) values,
%    format = 'uv'
%
%    If the user asks for just one variable returned, then we return a
%    2-vector that has uv(1:2), as in
%
%    N.B. There are two very closely related (u, v) formats.
%    These are (u, v) and (u', v'). The relationship between them is
%    uprime = u and vprime = 1.5 * v
%
%    The u', v' was an improvement made in the 1960s or so, and this
%    routine (by default) returns those values. The reason is because this
%    is mainly what people want. This function returns the (u'v') value
%    because, well, it is more modern and it is used in the LUV format.
%
%    X + Y + Z = 0 is returned as u = v = 0.
%
%    This function contains examples of usage inline. To access these, type
%    'edit xyz2uv.m' into the Command Window.
%
% Inputs:
%    xyz    - Vector. A 1x3 vector of the normalized CIE XYZ coordinates.
%    format - (Optional) String. This is how to distinguish between [u, v],
%             & [u', v']. The only case we care about is 'uv'. Default ''.
%
% Outputs:
%    u      - Numeric. The u-prime chromacity coordinate.
%    v      - Numeric. The v-prime chromacity coordinate.
%
% Optional key/value pairs:
%    None.
%
% References:
%    See (e.g.) Wyszecki and Stiles, 2cd, page 165.
%    <http://en.wikipedia.org/wiki/CIELUV>
%
% See Also:
%   cct, spd2cct, xyz2srgb, xyz2<TAB>
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    11/01/17  jnm  Comments & formatting
%    11/17/17  jnm  Formatting
%    07/15/19  JNM  Formatting update

% Examples:
%{
    wave = 400:10:700;
    d65 = ieReadSpectra('D65', wave);
    XYZ = ieXYZFromEnergy(d65', wave);
    [uP, vP] = xyz2uv(XYZ)
    u = uP;
    v = vP / 1.5;

    [u, v] = xyz2uv(XYZ, 'uv')
%}

if notDefined('xyz'), error('XYZ values required'); end
if notDefined('format'), format = ''; end
if size(xyz, 2) ~= 3, error('XYZ should be n x 3'); end

% Now compute uprime and vprime. These are the pre-cursors to ustar and
% vstar. The columns of xyz are X, Y and Z respectively.
B = (xyz(:, 1) + 15 * xyz(:, 2) + 3 * xyz(:, 3));

u = zeros(size(xyz, 1), 1);
v = zeros(size(u));

% Whenever B is valid, we set the u, v values to something legitimate. I am
% not sure what they should be when X + Y + Z is zero, as above. For now,
% we are leaving them as zero.
nz = B > 0;
u(nz) = 4 * xyz(nz, 1) ./ B(nz);
v(nz) = 9 * xyz(nz, 2) ./ B(nz);

% Check if the old 1960s (u, v), not (u', v') is being requested.
if isequal(format, 'uv'), v = v / 1.5; end

% If the user asked for just one ouptut, we combine u and v into a vector
% and return that
if nargout == 1
    u = [u, v];
end

end