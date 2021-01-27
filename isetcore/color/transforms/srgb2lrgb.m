function lrgb = srgb2lrgb(srgb)
% Transform srgb (nonlinear) to linear rgb (lrgb)
%
% Syntax:
%   lrgb = srgb2lrgb(srgb)
%
% Description:
%    sRGB is a standard color space that is designed around the properties
%    of display monitors (Sony Trinitron CRT). The RGB coordinates in this
%    space are nonlinearly related to XYZ. Prior to the nonlinear step,
%    however, there is a linear stage that we refer to as linear rgb
%    (lrgb). This routine converts the nonlinear (srgb) values into the
%    desired (lrgb) values.
%
%    The input range for srgb values is (0, 1); the range for the linear
%    values is (0, 1).
%
%    This function contains examples of usage inline. To access these, type
%    'edit srgb2lrgb.m' into the Command Window.
%
% Inputs:
%    srgb - The standard SRGB coordinates, prior to manipulation.
%
% Outputs:
%    lrgb - The linear RGB coordinates, prior to the sRGB scaling and
%           exponent
% Optional key/value pairs:
%    None.
%
% References:
%    <http://en.wikipedia.org/wiki/SRGB>
%    <http://www.w3.org/Graphics/Color/sRGB (techy)>
%
% See Also:
%   lrgb2srgb

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    11/01/17  jnm  Comments & formatting
%    11/17/17  jnm  Formatting
%    12/12/17  baw  Changed variable names to srgb/lrgb to respond to note
%    12/21/17  dhb  Make sure output dimensionality matches input, to
%                   preserve RGB or XW format from input to output.
%    07/15/19  JNM  Formatting update

% Examples:
%{
    %  The (1,1,1) RGB value maps to the XYZ of a D65 display with unit
    % luminance. This is XYZ = (0.9505, 1.0000, 1.0890)
    m = colorTransformMatrix('srgb2xyz');
    e = ones(1, 3) * m;
    e = round(e * 1000) / 1000
%}
%{
    srgb = [1, 1, 1] * 0.5;
    lrgb = srgb2lrgb(srgb)
    lrgb2srgb(lrgb)
%}

if max(srgb(:)) > 1
    warning('srgb2lrgb: srgb appears to be outside the (0,1) range');
end

% Allocate lrgb as same size as srgb. Need to to this so that format
% (RGB or XW) is preserved by this routine.
lrgb = zeros(size(srgb));

% Change to linear rgb values
big = srgb > 0.04045;
lrgb(~big) = srgb(~big) / 12.92;
lrgb(big) = ((srgb(big) + 0.055) / 1.055) .^ 2.4;

end
