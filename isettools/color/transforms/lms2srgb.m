function srgb = lms2srgb(lms)
% Convert LMS data to sRGB format for visualization
%
% Syntax:
%   srgb = lms2srgb(lms)
%
% Description:
%    Convert LMS data to sRGB format for visualization
%
% Inputs:
%    lms  - LMS Data - (what format?)
%
% Outputs:
%    srgb - standard Red-Green-Blue format
%
% See Also:
%    s_HumanColorBlind
%
% (c) ImagEval copyright 2012

% Examples:
%{
  scene    = sceneCreate; imgXYZ   = sceneGet(scene,'xyz');
  whiteXYZ = sceneGet(scene,'illuminant xyz');

  lms = xyz2lms(imgXYZ, 1, 'Brettel', whiteXYZ);  % Protan view
  imagesc(lms2srgb(lms))

  lms = xyz2lms(imgXYZ, 0);  % Normal view
  imagesc(lms2srgb(lms))
%}

srgb = xyz2srgb(imageLinearTransform(lms, colorTransformMatrix('lms2xyz')));

end