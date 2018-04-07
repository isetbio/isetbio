function img  = imgSweep(imSize, maxFreq)
% Create a sweep frequency image as a test pattern.
%
% Syntax:
%	img  = imgSweep([imSize], [maxFreq])
%
% Description:
%    The frequency increases across the columns; the contrast is high at
%    the top row and decreases down the rows. Used by sceneWindow.
%
%    There are examples contained in the code. To access, type 'edit
%    imgSweep.m' into the Command Window.
%
% Inputs:
%    imSize  - (Optional) The image size. Default 128.
%    maxFreq - (Optional) The maximum frequency. Default imSize/16 (8).
%
% Outputs:
%    img     - The created image
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Programming Note: We should be able to pass in the local frequency
%      as a function of image size, as well as the yContrast as a function
%      of image size. So, perhaps the function should be
%           img = imgSweep(imSize, [xFreq], [yContrast])
%      where length(xFreq) = length(yContrast) = imSize.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/29/18  jnm  Formatting

% Examples:
%{
    img  = imgSweep(256, 16);
    imagesc(img);
    colormap(gray);
    axis image

    img =  imgSweep;
    imagesc(img);
    colormap(gray);
    axis image
%}
if notDefined('imSize'), imSize = 128; end
if notDefined('maxFreq'), maxFreq = imSize / 16; end

% X positions in the image.
x = (1:imSize) / imSize;

% The change in frequency is slow at first, and then reaches a maximum
% frequency of one quarter the image size.
freq = (x .^ 2) * maxFreq;
xImage = sin(2 * pi * (freq .* x));
yContrast = (imSize:-1:1) / imSize;

img = yContrast' * xImage + 0.5;
img =  ieScale(img, 0, 255);

end