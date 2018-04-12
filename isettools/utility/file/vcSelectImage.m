function [fullName, imageType] = vcSelectImage(imageType, ext)
% Return the full path name and type of an ISET image.
%
% Syntax:
%   [fullName, imageType] = vcSelectImage([imageType], [ext])
%
% Description:
%    Interface to select an image from the isetbioDataPath image
%    directory. The imageType specifies the color type of the image
%    (monochrome, multispectral, rgb, hyperspectral).
%
% Inputs:
%    imageType - (Optional) The image color type. Options are monochrome,
%                multispectral, rgb, and hyperspectral. Default is ''
%    ext       - (Optional) The file extension type. Default is '*'
%
% Outputs:
%    fullName  - The full path and file name.
%    imageType - The image color type. Options are the subdirectory names
%                of isetbioDataPath/images, 'targets', 'multispectral',
%                'rgb', 'unispectral'.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: We should dir() the directory and verify that the string is
%      one of these.

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    11/29/17  jnm  Formatting, added note
%    01/17/18  dhb  Don't autorun the example
%    01/29/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip.  Requires user input.
    fullName = vcSelectImage()
    fullName = vcSelectImage('rgb')
    fullName = vcSelectImage('unispectral')
    fullName = vcSelectImage('multispectral')
%}

if notDefined('imageType'), imageType = ''; end
if notDefined('ext'), ext = '*'; end

if notDefined('imgDir')
    imgDir = fullfile(isetbioDataPath, 'images', imageType);
end

fullName = vcSelectDataFile(imgDir, 'r', ext);

if isempty(fullName)
    imageType = '';
    disp('User canceled');
elseif nargout == 2
    % Try to determine the image type, if this is requested.
    imageType = ieImageType(fullName);
end

end
