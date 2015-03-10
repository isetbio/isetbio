function [fullName, imageType] = vcSelectImage(imageType,imgDir,ext)
% Return the full path name and type of an ISET image.
%
%   [fullName, imageType] = vcSelectImage([imageType],[imgDir],[ext])
%
%  The FULLNAME contains the full path.  The imageType specifies the color
%  type of the image (monochrome, multispectral, rgb, hyperspectral).
%
%  If imgDir is not set, the routine assumes the file is in the ISET images
%  directory.
%
% Examples:
%  fullName = vcSelectImage('rgb')
%  fullName = vcSelectImage('rgb',pwd)
%  fullName = vcSelectImage('monochrome')
%  fullName = vcSelectImage('multispectral')
%
% Copyright ImagEval Consultants, LLC, 2003.

% Programming notes:
%  I am concerned whether the imageType is always determined correctly.

if notDefined('imageType'), imageType = ''; end
if notDefined('ext'), ext = '*'; end

curDir = pwd;
if notDefined('imgDir')
    imgDir = fullfile(isetRootPath,'data','images',imageType);
end

chdir(imgDir)
fullName = vcSelectDataFile('stayput','r',ext);

if isempty(fullName)
    imageType = '';
elseif nargout == 2
    % Try to determine the image type, if this is requested.
    imageType = ieImageType(fullName);
end

chdir(curDir)

end

