function [hc, rect] = hcimageCrop(img, rect, cPlane)
% Select region to crop from a hypercube image
%
% Syntax:
%   [hc, rect] = hcimageCrop(img, rect, cPlane)
%
% Description:
%    Select region to crop from a hypercube image
%
%    Examples are located within the code. To access the examples, type
%    'edit hcimageCrop.m' into the Command Window.
%
% Inputs:
%    img    - Hypercube input (required)
%    rect   - (Optional) If you know the rect, send it in. Default is to
%             show a compressed image in a window created by Matlab's
%             imcrop function, and get that to return the rect.
%    cPlane - (Optional) Which plane to use for cropping. Default 1
%
% Outputs:
%    hc     - Cropped hc image
%    rect   - rectangle used for cropping
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    hcimage, hcimageRotateClip


% History:
%    xx/xx/xx       (c) Imageval
%    12/06/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip.  Requires user input.
    %
    % To get this to work, you need to understand now
    % function imcrop works.
    mcc = macbethChartCreate;
    [hc, rect] = hcimageCrop(mcc.data, [], 1);
%}

if notDefined('img'), error('hyper cube image required'); end
if notDefined('cPlane'), cPlane = 1; end

if notDefined('rect')
    %Use a square root to avoid problems with bright, saturating pixels.
    tmp = sqrt(double(img(:, :, cPlane)));
    tmp = tmp / max(tmp(:));
    f = vcNewGraphWin;
    imagesc(tmp);
    colormap(hot);
    [d, rect] = imcrop;
    close(f)
else
    d = imcrop(double(img(:, :, cPlane)), rect);
end

% Create space for cropped image
[r, c] = size(d);
w = size(img, 3);
if isa(img, 'uint16')
    hc = zeros(r, c, w, 'uint16');
else
    hc = zeros(r, c, w, 'double');
end

% Crop each plane
h = waitbar(0, 'Cropping');
for ii = 1:w
    waitbar(ii / w, h);
    hc(:, :, ii) = imcrop(img(:, :, ii), rect);
end
close(h);

end