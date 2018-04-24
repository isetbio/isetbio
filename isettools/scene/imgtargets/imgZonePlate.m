function img = imgZonePlate(sz, amp, ph)
% Make a zone plate image
%
% Syntax:
%   img = ieZonePlate([SZ], [AMP], [PHASE])
%
% Description:
%    Create a zone plate image, which follows the function:
%
%       AMP * cos( r ^ 2 + PHASE) + 1
%
%    There are examples contained in the code. To access, type 'edit
%    imgZonePlate.m' into the Command Window.
%
% Inputs:
%    sz  - (Optional) The image size. Default 256.
%    amp - (Optional) The image amplitude. Default is 1.
%    ph  - (Optional) The image phase. Default 0.
%
% Outputs:
%    img - The created image
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/29/18  jnm  Formatting

% Examples:
%{
    img = imgZonePlate(256);
    imagesc(img);
    axis image;
    colormap(gray(256)); 
%}
%{
    img = imgZonePlate([384, 384], 255);
    imagesc(img);
    colormap(gray(256));
    axis image
%}
%{
    img = imgZonePlate([384, 384], 255, pi / 2);
    imagesc(img);
    colormap(gray(256));
    axis image
%}

if notDefined('sz'),  sz = [256 256]; end
if notDefined('amp'), amp = 1; end
if notDefined('ph'),  ph = 0;   end

if (length(sz) == 1),  sz = [sz, sz]; end

mxsz = max(sz(1), sz(2));

% Almost ordinary, except we adjus to a minimum value of 0
img = amp * cos( (pi / mxsz) * imgRadialRamp(sz, 2) + ph ) + 1;
% img = ieScale(img, 0, 1);

end
