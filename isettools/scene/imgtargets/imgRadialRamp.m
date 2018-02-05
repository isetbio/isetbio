function img = imgRadialRamp(sz, expt, origin)
% Make a radial ramp function
%
% Syntax:
%	img = imgRadialRamp(sz, expt, origin)
%
% Description:
%    The radial ramp can be either linear (default) or raised to power EXPT
%    (default = 1).  The origin is normally in the center of the image
%    (default = (size + 1) / 2, [1 1] = upper left), but its location can
%    be specified (ORIGIN).  The default size sz = (row, col) is 256, 256.
%
%    There are examples contained in the code. To access, type 'edit
%    imgRadialRamp.m' into the Command Window.
%
% Inputs:
%    sz     - (Optional) The [row, col] size of the image.
%             Default is [128, 128].
%    expt   - (Optional) The exponential. Default is 1.
%    origin - (Optional) The specified location of the origin.
%             Default is [1, 1] (upper left)
%
% Outputs:
%    img    - The created radial ramp function image
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    02/02/18  jnm  Formatting

% Examples:
%{
    img = imgRadialRamp([384, 384]);
    imagesc(img);
    axis image
%}
%{
    img = imgRadialRamp([128, 128], 2, [1, 1]);
    imagesc(img);
    axis image
%}

if notDefined('sz'), sz = [128, 128]; end
if notDefined('expt'), expt = 1;  end
if notDefined('origin'), origin = (sz+1)/2; end

[xramp, yramp] = meshgrid( (1:sz(2)) - origin(2), (1:sz(1)) - origin(1) );

img = (xramp .^ 2 + yramp .^ 2) .^ (expt / 2);

end