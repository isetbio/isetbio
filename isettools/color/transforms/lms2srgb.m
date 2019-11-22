function srgb = lms2srgb(lms)
% Convert LMS data to sRGB format for visualization
%
% Syntax:
%   srgb = lms2srgb(lms)
%
% Description:
%    Convert LMS data to sRGB format for visualization
%
%    This function contains examples of usage inline. To access these, type
%    'edit lms2srgb.m' into the Command Window.
%
% Inputs:
%    lms  - Matrix. LMS Data - (what format?)
%
% Outputs:
%    srgb - Matrix. A standard Red-Green-Blue format image matrix.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   s_HumanColorBlind
%

% History:
%    XX/XX/12       (c) ImagEval copyright 2012
%    07/16/19  JNM  Formatting

% Examples:
%{
    scene = sceneCreate;
    imgXYZ = sceneGet(scene,'xyz');
    whiteXYZ = sceneGet(scene,'illuminant xyz');

    lms = xyz2lms(imgXYZ, 1, 'Brettel', whiteXYZ);  % Protan view
    imagesc(lms2srgb(lms))

    lms = xyz2lms(imgXYZ, 0);  % Normal view
    imagesc(lms2srgb(lms))
%}

srgb = xyz2srgb(imageLinearTransform(lms, ...
    colorTransformMatrix('lms2xyz')));

end