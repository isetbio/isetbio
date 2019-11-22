function srgb = lrgb2srgb(lrgb)
% Convert linear sRGB values to proper sRGB values
%
% Syntax:
%   srgb = lrgb2srgb(lrgb)
%
% Description:
%    This routine implements the nonlinear step for converting linear rgb
%    (lrgb) into the frame buffer representations in the srgb
%    representation. The rgb data can be in either RGB or XW format.
%
%    The inputs are linear rgb values, and the returned values are
%    nonlinear framebuffer values. They are in the range [0, 1]
%
%    The gamValue used in the srgb formula combines with a linear regime
%    that makes an overall approximation of the display gamma as 2.2
%
%    This function contains examples of usage inline. To access these, type
%    'edit lrgb2srgb.m' into the Command Window.
%
% Inputs:
%    lrgb - Vector. The linear RGB coordinates, prior to manipulation.
%
% Outputs:
%    srgb - Vector. The standard RGB coordinates, after manipulation.
%
% Optional key/value pairs:
%    None.
%
% References:
%    Previous web references have been deleted from this file because the
%    links are now dead. The proper link is to the Wikipedia srgb page.
%
%    <https://en.wikipedia.org/wiki/SRGB>
%
% See Also:
%   xyz2srgb, srgb2lrgb
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC.
%    11/01/17  jnm  Comments & formatting
%    12/22/17  dhb  Don't use same variable name for input and output
%    07/16/19  JNM  Formatting update, add example

% Examples:
%{
    lrgb = [1, 1, 1] / 2;
    srgb = lrgb2srgb(lrgb)
%}

if max(lrgb(:)) > 1 || min(lrgb(:)) < 0
    error('Linear rgb values must be between 0 and 1');
end

% These are framebuffer values, but they live in the [0, 1] space. The
% transformation, which has a linear and power part, is intended to
% approximate a gamma of 2.2 as a whole.
srgb = zeros(size(lrgb));
big = lrgb > 0.0031308;
srgb(~big) = lrgb(~big) * 12.92;
srgb(big) = 1.055 * lrgb(big) .^ (1 / 2.4) - 0.055;

end
