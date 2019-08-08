function [XW, r, c, w] = RGB2XWFormat(imRGB)
% Transform an RGB form matrix into an XW (space-wavelength) matrix
%
% Syntax:
%   [XW, r, c, w] = RGB2XWFormat(imRGB)
%
% Description:
%    This  routine converts from RGB format to XW format. The row, column
%    and w (number of color bands) of the imRGB are also returned, if
%    requested.
%
%    We say matrices in (r, c, w) format are in RGB format. The dimension,
%    w, represents the number of data color bands. When w = 3, the data
%    are an RGB image. But w can be almost anything (e.g., 31 wavelength
%    samples from 400:10:700). We will use this format frequently for
%    spectral data.
%
%    The RGB format is useful for imaging. When w = 3, you can use
%    conventional image() routines. When w > 3, use imageSPD.
%
%    The XW (space-wavelength) format is useful for computation. In this
%    format, for example, XW * spectralFunction yields a spectral response.
%
%    The inverse routine is XW2RGBFormat
%
%    This function contains examples of usage. To access, type 'edit
%    RGB2XWFormat.m' into the Command Window.
%
% Inputs:
%    imRGB - Matrix. The provided RGB formatted matrix.
%
% Outputs:
%    XW    - Matrix. The space-wavelength formatted matrix
%    r     - Numeric. Row data for the space-wavelength matrix (# rows).
%    c     - Numeric. Column data for the space-wavelength matrix (# cols).
%    w     - Numeric. The number of data color bands.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   imageSPD, imagescRGB, XW2RGBFormat
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/03/19  JNM  Formatting update

% Examples:
%{
    xwSRGBs = [[188 188 188]' [124 218 89]' [255 149 203]' ...
       [255 3 203]']' / 255;
    rgbSRGBs = XW2RGBFormat(xwSRGBs, 2, 2);
    [xwCheck, r, w, c] = RGB2XWFormat(rgbSRGBs)
%}
s = size(imRGB);

% If the data are in a matrix, then assume only one wavelength dimension,
% (row, col, 1).
if length(s) < 3, s(3) = 1; end
XW = reshape(imRGB, s(1) * s(2), s(3));

r = s(1);
c = s(2);
w = s(3);

end