function imIndices = rgcGetIndicesOfSubImage(imgPixSize,subImgCenter,subImgSize,coneSpacing)
% Computes the indices corresponding to the subImage
% 
%   imIndices =
%   rgcGetIndicesOfSubImage(imgSize,subImgCenter,subImgSize,coneSpacing);
%
% (c) Stanford Synapse Team 2010

imMin = subImgCenter-(subImgSize-1)/2;
imMax = subImgCenter+(subImgSize-1)/2;

ImgCenter = (imgPixSize+1)*coneSpacing/2;

% position in the space centered at image(1,1) (um)
posMin = imMin+ImgCenter;
posMax = imMax+ImgCenter;

% same positions in pixels
posMin = floor(posMin/coneSpacing);
posMax = ceil(posMax/coneSpacing);

imIndices{1} = posMin(1):posMax(1);
imIndices{2} = posMin(2):posMax(2);