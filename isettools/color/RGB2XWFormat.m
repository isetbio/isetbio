function [XW, r, c, w] = RGB2XWFormat(imRGB)
% Transform an RGB form matrix into an XW (space-wavelength) matrix
%
% Syntax:
%   [XW, r, c, w] = RGB2XWFormat(imRGB)
%
% Description:
%    This  routine converts from RGB format to XW format.  The row and
%    column of the imRGB are also returned, if requested.
%
%    We say matrices in (r, c, w) format are in RGB format.  The dimension,
%    w, represents the number of data color bands.  When w=3, the data are
%    an RGB image. But w can be almost anything (e.g., 31 wavelength
%    samples from 400:10:700).  We will use this format frequently for
%    spectral data.
%
%    The RGB format is useful for imaging.  When w = 3, you can use
%    conventional image() routines.  When w > 3, use imageSPD.
%
%    The XW (space-wavelength) format is useful for computation.  In this
%    format, for example, XW*spectralFunction yields a spectral response.
%
%    The inverse routine is XW2RGBFormat
%
% Inputs:
%    imRGB - The provided RGB formatted matrix
%
% Outputs:
%    XW    - Space-wavelength formatted matrix
%    r     - Row data for the space-wavelength matrix
%    c     - Column data for the space-wavelength matrix
%    w     - The number of data color bands
%
% See Also:
%    imageSPD, imagescRGB, XW2RGBFormat
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting

% Examples:
%{
   ptbSRGBs = [[188 188 188]' [124 218 89]' [255 149 203]' [255 3 203]'];

   % The ISET form takes the frame buffer values in the [0,1] regime
   isetSRGBs = ptbSRGBs/255;
   isetSRGBs = XW2RGBFormat(isetSRGBs',4,1);
   isetXYZ   = srgb2xyz(isetSRGBs);
   isetXYZs  = RGB2XWFormat(isetXYZ)';
%}
s = size(imRGB);

% If the data are in a matrix, then assume only one wavelength dimension, 
% (row, col, 1).
if length(s) < 3
    s(3) = 1;
end

XW = reshape(imRGB, s(1)*s(2), s(3));

r = s(1);
c = s(2);
w = s(3);

end