function [hc, rect] = hcimageCrop(img, rect, cPlane)
% Select region to crop from a hypercube image
%
% Syntax:
%   [hc, rect] = hcimageCrop(img, rect, cPlane)
%
% Description:
%    Select region to crop from a hypercube image
%
% Inputs:
%    img    - Hypercube input (required)
%    rect   - (Optional) If you know the rect, send it in. Default is to
%             use cPlane in a square root to avoid problems with bright,
%             saturating pixels when creating the rect.
%    cPlane - (Optional) Which plane to use for cropping. Default 1
%
% Outputs:
%    hc     - Cropped hc image
%    rect   - rectangle used for cropping
%
% Notes:
%    * [Note: JNM - Example does not work!]

% History:
%    xx/xx/xx       (c) Imageval
%    12/06/17  jnm  Formatting

% Examples:
%{
    img = macbethChartCreate;
    [hc, rect] = hcimageCrop(img, [], 1)
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
for ii=1:w
    waitbar(ii / w, h);
    hc(:, :, ii) = imcrop(img(:, :, ii), rect);
end
close(h);

end