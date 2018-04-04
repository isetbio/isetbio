function img  = imgMackay(radialFrequency, imSize)
% Create a MacKay chart spatial pattern. 
%
% Syntax:
%	img  = imgMackay([radialFrequency], [imSize])
%
% Description:
%    The Mackay chart has lines at many angles and increases in spatial
%    frequency from periphery to center. This routine is called in creating
%    the Mackay scene (sceneWindw).
%
%    There are examples contained in the code. To access these, simply type
%    'edit imgMackay.m' into the Command Window.
%  
% Inputs:
%    radialFrequency - (Optional) The radial frequency of the image.
%                      Default is 8.
%    imSize          - (Optional) The image size. Default is 128.
%
% Outputs:
%    img             - The created image.
%
% Optional key/value pairs:
%	 None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    02/02/18  jnm  Formatting

% Examples:
%{
    img =  imgMackay;
    imagesc(img);
    colormap(gray);
    axis image
%}
%{
    img  = imgMackay(12, 256);
    imagesc(img);
    colormap(gray);
    axis image
%}

if notDefined('radialFrequency'), radialFrequency = 8; end
if notDefined('imSize'), imSize = 128; end

mx = round(imSize / 2);
mn = -(mx - 1);
[x, y] = meshgrid(mn:mx, mn:mx);
l = (x == 0);
x(l) = eps;

img = cos(atan(y ./ x) * 2 * radialFrequency);
img =  ieScale(img, 0, 255);

end