function lSize = convertPixelsInMicrons(subPatchesSize,totalSize,fov,scene,oi)
% Converts the size in number of pixels, to a size in um 
% (that could be on the retina, or on the rgc)
%
%   lSize = convertPixelsInMicrons(subPatchesSize,totalSize,fov,scene,oi);
%
% subPatchesSize: Size to be converted
% totalSize: Size of the whole image
% fov : field of view for the whole image.
% scene : ISET scene
% oi : ISET oi
%
% (c) Stanford Synapse Team - 2010

warning('Should consider refine this function. Not called currently');

% number of microns per degree
micronPerDegree = rgcGetMicronPerDegree(scene,oi);

% (proportion of the image) * (size of the image in microns)
lSize = (subPatchesSize/totalSize)*(fov*micronPerDegree);