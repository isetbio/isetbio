function imRGB = XW2RGBFormat(imXW, row, col)
% Convert XW format data to RGB format
%
% Syntax:
%   imRGB = XW2RGBFormat(imXW, row, col)
%
% Description:
%    This  routine converts from XW format to RGB format. The row and
%    column of the imXW are required input arguments.
%
%    We say matrices in (r, c, w) format are in RGB format. The dimension,
%    w, represents the number of data color bands. When w = 3, the data are
%    an RGB image. But w can be almost anything (e.g., 31 wavelength
%    samples from 400:10:700). We use this format frequently for spectral
%    data.
%
%    The RGB format is useful for imaging. When w = 3, you can use
%    conventional image() routines. When w > 3, use imageSPD.
%
%    The XW (space-wavelength) format is useful for computation. In this
%    format, for example, XW * spectralFunction yields a spectral response.
%
%    The inverse routine is RGB2XWFormat.
%
%    This function contains examples of usage inline. To access these, type
%    'edit XW2RGBFormat.m' into the Command Window.
%
% Inputs:
%    imXW  - Matrix. The space-wavelength formatted data.
%    row   - Numeric. The imXW row data (number of rows).
%    col   - Numeric. The imXW column data (number of columns).
%
% Outputs:
%    imRGB - Matrix. The RGB format data.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   imageSPD, imagescRGB, RGB2XWFormat
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC
%    10/27/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/03/19  JNM  Formatting update

% Examples:
%{
   xwSRGBs = [[188 188 188]' [124 218 89]' [255 149 203]' ...
       [255 3 203]']' / 255;
   rgbSRGBs = XW2RGBFormat(xwSRGBs, 2, 2);
%}

if notDefined('imXW'), error('No image data.'); end
if notDefined('row'), error('No row size.'); end
if notDefined('col'), error('No col size.'); end

x = size(imXW, 1);
w = size(imXW, 2);

if row * col ~= x, error('XW2RGBFormat:  Bad row, col values'); end

imRGB = reshape(imXW, row, col, w);

end