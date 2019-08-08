function [srgb, lrgb, maxY] = xyz2srgb(xyz)
% Convert CIE XYZ to sRGB color space
%
% Syntax:
%   [srgb, lrgb, maxY] = xyz2srgb(xyz)
%
% Description
%    The CIE XYZ values are in an RGB Format image. They are converted to
%    sRGB values. The user can also get linear RGB values as well as the Y
%    value (maxY) that is used to scale the XYZ image so that it is within
%    the [0, 1] range as required by the sRGB standard.
%
%    The sRGB color space is a display-oriented representation that matches
%    a Sony Trinitron. The monitor white point is assumed to be D65. The
%    white point chromaticity are (.3127, .3290), and for an sRGB display
%    (1, 1, 1) is assumed to map to XYZ = (0.9504 0.9999 1.0891).
%    The RGB primaries of an srgb display have xy coordinates of
%    xy = [.64, .3; .33, .6; .15, .06]
%
%    The overall gamma of an sRGB display is about 2.2, but this is because
%    at low levels the value is linear and at high levels the gamma is 2.4.
%    See the wikipedia page for a discussion.
%
%    sRGB values run from [0 1]. At Imageval this assumption changed from
%    the range [0 255] on July 2010. This was based on the wikipedia entry
%    and discussions with Brainard. Prior calculations of delta E are not
%    changed by this scale factor.
%
%    The linear srgb values (lRGB) can also be returned. These are the
%    values of the linear phosphor intensities, without any gamma or
%    clipping applied. lRGB values nominally run from [0, 1], but we allow
%    them to be returned  outside of this range.
%
%    This function contains examples of usage inline. To access these, type
%    'edit xyz2srgb.m' into the Command Window.
%
% Inputs:
%    xyz  - Matrix. The XYZ values, in isetbio RGB image format.
%
% Outputs:
%    srgb - Matrix. The standard Red-Green-Blue values.
%    lrgb - Matrix. The linear Red-Green-Blue values.
%    maxY - Numeric. The Y value, which can be used to scale the xyz image
%           to the sRGB standard.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    This xyz -> sRGB matrix is supposed to work for XYZ values scaled so
%    that the maximum Y value is around 1. In the Wikipedia page, they say:
%       if you start with XYZ values going to 100 or so, divide them by 100
%       first, or apply the matrix and then scale by a constant factor to
%       the [0, 1] range).
%    They add:
%       display white represented as (1, 1, 1) [RGB]; the corresponding
%       original XYZ values are such that white is D65 with unit luminance
%       (X, Y, Z = 0.9505, 1.0000, 1.0890).
%
% References:
%    Modern: <http://en.wikipedia.org/wiki/SRGB>
%    Original: <http://www.w3.org/Graphics/Color/sRGB>
%
% See Also:
%   srgb2xzy, lrgb2srgb, colorTransformMatrix, mageLinearTransform.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    11/01/17  jnm  Comments & formatting
%    11/17/17  jnm  Formatting
%    07/15/19  JNM  Formatting update

% Examples:
%{
    inputSRGBs = [[188 188 188]' [124 218 89]' [255 149 203]' ...
        [255 3 203]'] / 255;
    isetSRGBs = XW2RGBFormat(inputSRGBs', 4, 1)
    isetXYZ = srgb2xyz(isetSRGBs);
    againSRGBs = xyz2srgb(isetXYZ)
%}

% The matrix converts (R, G, B) * matrix. This is the transpose of the
% Wikipedia page.
matrix = colorTransformMatrix('xyz2srgb');

% Notice that (1, 1, 1) maps into D65 with unit luminance (Y)
% matrix = colorTransformMatrix('srgb2xyz');
% ones(1, 3) * matrix

% The linear transform is built on the assumption that the maximum
% luminance is 1.0. If the inputs are all within [0, 1], I suppose we
% should leave the data alone. If the maximum XYZ value is outside the
% range, we need to scale. We return the true maximum luminance in the
% event the user wants to invert, later.
Y = xyz(:, :, 2);
maxY = max(Y(:));
if maxY > 1
    xyz = xyz / maxY;
else
    maxY = 1;
end
if min(xyz(:)) < 0
    fprintf('Warning:  Clipping negative values in XYZ %f\n', min(xyz(:)));
    xyz = ieClip(xyz, 0, 1);
end
lrgb = imageLinearTransform(xyz, matrix);

% The sRGB values must be clipped to 0, 1 range.
% The linear values may be outside the range. This is also described on
% the Wikipedia page.
srgb = lrgb2srgb(ieClip(lrgb, 0, 1));

end