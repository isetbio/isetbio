function img  = imgRamp(imSize)
% Create a set of intensity ramps as test spatial pattern. 
%
% Syntax:
%	img  = imgRamp([imSize])
%
% Description:
%    Ramp patterns are useful test patterns (sceneWindow) for evaluating
%    contouring caused by poor analog to digital conversion, and sometimes
%    for evaluating problems with demosaic'ing routines. This routine creates
%    a steep intensity ramp at the top row of the image, and a decreasing
%    intensity ramp as we measure down the rows.
%
%    There are examples contained in the code. To access, type 'edit
%    imgRamp.m' into the Command Window.
%
% Inputs:
%    imSize - (Optional) The image size. Default 128.
%
% Outputs:
%    img    - The created image
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    02/02/18  jnm  Formatting

% Examples:
%{
    img = imgRamp(128);
    imagesc(img);
    colormap(gray);
    axis image
    mesh(img)
%}
%{
    img = imgRamp;
    imagesc(img);
    colormap(gray);
    axis image
    imagesc(img);
    colormap(gray(256));
    truesize;
%}

if notDefined('imSize'), imSize = 128; end

% X positions in the image.
mx = round(imSize / 2);
mn = -(mx - 1);
xImage = mn:mx;

yContrast = ((imSize:-1:1) / imSize);
img = (yContrast' * xImage) + 0.5;
img =  ieScale(img, 0, 255);

end