function img  = imgRamp(imSize)
% Create a set of intensity ramps as test spatial pattern.  
%
%   img  = imgRamp(imSize)
%
% Ramp patterns are useful test patterns (sceneWindow) for evaluating
% contouring caused by poor analog to digital conversion, and sometimes
% for evaluating problems with demosaic'ing routines.  This routine creates
% a steep intensity ramp at the top row of the image, and a decreasing
% intensity ramp as we measure down the rows.
%
% Examples:
%   figure; 
%   img  = imgRamp(128); imagesc(img); colormap(gray); axis image
%   mesh(img)
%   img =  imgRamp; imagesc(img); colormap(gray); axis image
%   imagesc(img); colormap(gray(256)); truesize;
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('imSize'), imSize = 128; end

% X positions in the image.
mx = round(imSize/2);
mn = -(mx-1);
xImage = mn:mx;

yContrast = ((imSize:-1:1)/imSize);
img = (yContrast'*xImage) + 0.5;
img =  ieScale(img,0,255);

end