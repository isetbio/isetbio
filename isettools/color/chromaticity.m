function xy = chromaticity(XYZ)
% Compute CIE chromaticity (xy) coordinates from XYZ data
%
% Syntax:
%   xy = chromaticity(XYZ)
%
% Description:
%    The input (XYZ) data can be in XW (space-wavelength) or RGB format.
%    In XW format, we expect N rows corresponding to spatial positions and
%    three columns containing in X, Y and Z. The chromaticity coordinates
%    (x, y) are returned in the columns of an Nx2 matrix.
%
%    If the data are in RGB format, the three planes should be (X, Y, Z)
%    images. The returned data are in as a two dimensional image format,
%    with each spatial position containing the corresponding (x, y) value.
%
%    There is nothing sacred about using this for XYZ -> xy. You can pass
%    any tristimulus coordinates and get the first two corresponding
%    chromaticity coordinates. So, for example RGB -> rg.
%
%    This function contains examples of usage inline. To access these, type
%    'edit chromaticity.m' into the Command Window.
%
% Inputs:
%    XYZ - Matrix. The CIE XYZ color space data. In XW or RGB format.
%
% Outputs:
%    xy  - Matrix. The CIE Chromaticity coordinates. The matrix is the
%          format of the number of coordinates by 2, for the (x, y) format.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/30/17  jnm  Comments & formatting
%    07/11/19  JNM  Formatting update

% Examples:
%{
    patchSize = 1;
    macbethChart = sceneCreate('macbeth', patchSize);
    p = sceneGet(macbethChart, 'photons');
    wave = sceneGet(macbethChart, 'wave');
    e = Quanta2Energy(wave, p);
    XYZ = ieXYZFromEnergy(e, wave);
    chromaticity(XYZ)

    XYZ = RGB2XWFormat(XYZ);
    chromaticity(XYZ)
%}

if ismatrix(XYZ)
    if size(XYZ, 2) ~= 3
        error('The XW input data should be in the columns of a Nx3 matrix')
    end

    ncols = size(XYZ, 1);
    xy = zeros(ncols, 2);

    s = sum(XYZ, 2);
    p = find(s ~= 0);
    xy(p, 1) = XYZ(p, 1) ./ s(p);
    xy(p, 2) = XYZ(p, 2) ./ s(p);

elseif ndims(XYZ) == 3
    [r, c, ~] = size(XYZ);
    xy = zeros(r, c, 2);

    s = XYZ(:, :, 1) + XYZ(:, :, 2) + XYZ(:, :, 3);
    xy(:, :, 1) = XYZ(:, :, 1) ./ s;
    xy(:, :, 2) = XYZ(:, :, 2) ./ s;
else
    error('Data must be either Nx3 or NxMx3 with XYZ in image planes');
end

end