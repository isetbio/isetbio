function rgb = srgb2lrgb(rgb)
% Convert srgb (nonlinear) to linear rgb (lrgb) the linear precursor to srgb   
%
% Syntax:
%   lrgb = srgb2lrgb(srgb)
%
% Description:
%    sRGB is a standard color space that is designed around the properties
%    of display monitors (Sony Trinitron CRT). The RGB coordinates in this
%    space are nonlinearly related to XYZ.  Prior to the nonlinear step,
%    however, there is a linear stage that we refer to as linear rgb
%    (lrgb). This routine converts the nonlinear (srgb) values into the
%    desired (lrgb) values.
%
%    The input range for srgb values is (0,1); the range for the linear
%    values is (0,1).
%
% Inputs:
%    rgb - The standard Red-Green-Blue coordinates, prior to manipulation.
%
% Outputs:
%    rgb - The linear Red-Green-Blue coordinates, after sRGB manipulation.
%
% Notes:
% * Imageval scaling changes as of July, 2010, as per discussions with
%   Brainard.
% * I am concerned that the input and output variables share a name
%
% References:
%    <http://en.wikipedia.org/wiki/SRGB>
%    <http://www.w3.org/Graphics/Color/sRGB (techy)>
%
% See Also:
%    lrgb2srgb.m, s_SRGB
%
% Copyright ImagEval Consultants, LLC, 2003.

% Examples:
%{
   %  The (1,1,1) RGB value maps to the XYZ of a D65 display with unit
   % luminance. This is XYZ = (0.9505, 1.0000, 1.0890)
   m = colorTransformMatrix('srgb2xyz'); 
   e = ones(1,3)*m;  e = round(e*1000)/1000
%}

if max(rgb(:)) > 1
    warning('srgb2lrgb: srgb appears to be outside the (0,1) range');
end

% Change to linear rgb values
big = (rgb > 0.04045);
rgb(~big) = rgb(~big)/12.92;
rgb(big) = ((rgb(big)+0.055)/1.055).^ (2.4);

end
