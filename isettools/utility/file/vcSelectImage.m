function [fullName, imageType] = vcSelectImage(imageType, imgDir, ext)
% Return the full path name and type of an ISET image.
%
% Syntax:
%   [fullName, imageType] = vcSelectImage([imageType], [imgDir], [ext])
%
% Description:
%    The fullName contains the full path. The imageType specifies the color
%    type of the image (monochrome, multispectral, rgb, hyperspectral).
%
%    If imgDir is not set, the routine assumes the file is in the ISET
%    images directory.
%
% Inputs:
%    imageType - (Optional) The image color type. Options are monochrome,
%                multispectral, rgb, and hyperspectral. Default is ''
%    imgDir    - (Optional) Image directory. Default is ISET images
%    directory
%    ext       - The file extension type. Default is '*'
%
% Outputs:
%    fullName  - The full path and file name.
%    imageType - The image color type. Options are monochrome,
%                multispectral, rgb, and hyperspectral.
%
% Notes:
%    * [Note: XXX - Programming notes: I am concerned whether the imageType
%      is always determined correctly.]

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    11/29/17  jnm  Formatting, added note

% Examples:
%{
    fullName = vcSelectImage()
    fullName = vcSelectImage('rgb')
    fullName = vcSelectImage('rgb', pwd)
    fullName = vcSelectImage('monochrome')
    fullName = vcSelectImage('multispectral')
%}

if notDefined('imageType'), imageType = ''; end
if notDefined('ext'), ext = '*'; end

curDir = pwd;
if notDefined('imgDir')
    imgDir = fullfile(isetbioDataPath, 'images', imageType);
end

chdir(imgDir)
fullName = vcSelectDataFile('stayput', 'r', ext);

if isempty(fullName)
    imageType = '';
    disp('User canceled');
elseif nargout == 2
    % Try to determine the image type, if this is requested.
    imageType = ieImageType(fullName);
end

chdir(curDir)

end
