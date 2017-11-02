function xyz = lms2xyz(lms)
% Convert stockman lms to xyz 10 deg
%
% Syntax:
%   xyz = lms2xyz(lms)
%
% Description:
%    Convert stockman lms to xyz 10 degree
%
% Inputs:
%    lms - Probably in RGB image format. *Should check.*
%
% Outputs:
%    xyz - XYZ Formatted image data
%
% Notes:
%    * [Note: XXX - Why is the Input variable format in question and not
%      already known?]
%    * [Note: JNM - Why is the XW format below noted as not thoroughly
%      debugged? Is there something we can do to address this?]
%
% See Also:
%    xyz2lms
%

% History:
%    xx/xx/12       (c) ImagEval
%    11/01/17  jnm  Comments, formatting, and adding actual example.

% Examples:
%{
   scene = sceneCreate('reflectance chart');
   vcAddAndSelectObject(scene); sceneWindow
   imgLMS = sceneGet(scene, 'lms');
   
   imgXYZ = lms2xyz(imgLMS);
   vcNewGraphWin;
   imagescRGB(lms2srgb(imgXYZ));
%}

if notDefined('lms'), error('lms required'); end

if ndims(lms) == 3
    % RGBW format
    xyz = imageLinearTransform(lms, colorTransformMatrix('lms2xyz'));
elseif ismatrix(lms)
    % XW format - Not debugged thoroughly
    if size(lms, 1) == 3 && size(lms, 2) ~= 3
        xyz = lms' * colorTransformMatrix('lms2xyz');
    elseif size(lms, 1) ~= 3 && size(lms, 2) == 3
        xyz = lms * colorTransformMatrix('lms2xyz');
    else
        error('Ambiguous lms shape');
    end
end

end
