function [h, rgbim] = imagescRGB(rgbim, varargin)
% Display a scaled RGB format image. 
%
% Syntax:
%   [h, rgbim] = imagescRGB(rgbim, [gamma]);
%   [h, rgbim] = imagescRGB(rgbim, row, col, [gamma])
%
% Description:
%    Prior to display negative values are clipped, and the clipped data are
%    scaled to a maximum of 1.
%
%	 If the exponent gamma is included, then rgbim .^ gamma are displayed;
% 
%    The routine accepts data in XW and RGB format. 
%        In XW format use:  imagescRGB(img, row, col, [gamma])
%        In RGB format use: imagescRGB(img, [gamma])
%
%    Examples are located within the code. To access the examples, type
%    'edit imagescRGB.m' into the Command Window.
%
% Inputs:
%    rgbim - The RGB Image
%    varargin - The other potential variables, depending on the format of
%               the incoming data. For an image in RGB format, there is the
%               optional variable gamma. For an XW image, the variables row
%               and col for rows and columns are required, with gamma
%               remaining an optional variable.
%           row   - The number of rows in the XW image
%           col   - The number of columns in the XW image
%           gamma - (Optional) The Luminance of the image, both formats.
%
% Outputs:
%    h        - The image handle
%    rgbim    - The image data
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: I am concerned about the ordering of the ^ gamma and the scale
%      operations. Perhaps scaling should be first, and then the gamma. As
%      things stand, we apply gamma and then scale. That applies here and
%      in other routines.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/08/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    foo = load('trees');
    [r, c] = size(foo.X);
    for ii = 1:3, rgb(:, :, ii) = reshape(foo.map(foo.X, ii), r, c); end

    rgbScaled = imagescRGB(rgb);
    rgbScaled = imagescRGB(rgb, 0.3);

    rgbXW = RGB2XWFormat(rgb);
    rgbScaled = imagescRGB(rgbXW, r, c, 0.3);
%}

% This is a theory of display. I am not sure I should be clipping before
% scaling. But over the years, that has seemed better.
rgbim = ieClip(rgbim, 0, []);
s = max(rgbim(:));
if s ~= 0, rgbim = rgbim / max(rgbim(:)); end

if ismatrix(rgbim)
    if nargin < 3
        error('XW input requires row and col arguments.');
    else
        row = varargin{1};
        col = varargin{2};
    end
    
    rgbim = XW2RGBFormat(rgbim, row, col);
    if nargin > 3
        gamma = varargin{3};
        rgbim = rgbim .^ gamma;
    end
    
elseif ndims(rgbim) == 3
    % row = size(rgbim, 1); 
    % col = size(rgbim, 2);
    if nargin > 1
        gamma = varargin{1};
        rgbim = rgbim .^ gamma;
    end
else 
    error('Bad image input');
end

% Eliminated imshow and replaced with this so it would work on a Jupyter
% hub site.
h = image(rgbim);
axis image;
axis off

end