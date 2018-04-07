function scene = sceneAdjustLuminance(scene, meanL)
% Scale scene mean luminance to meanL
%
% Syntax:
%   scene = sceneAdjustLuminance(scene, meanL)
%
% Description:
%    The photon level in the scene structure is multiplied so that the
%    mean luminance level is meanL, rather than the current level. The
%    illuminant is also scaled to preserve the reflectance.
%
%    There are examples in the code. Type 'edit sceneAdjustLuminance' into
%    the Command Window to access.
%
% Inputs:
%    scene - A scene structure
%    meanL - The meanLuminance to scale the scene structure to.
%
% Outputs:
%    scene - The modified scene strucuture
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO - Some scenes go to very long wavelengths. That slows the
%      calculation. Never let the calculation go beyond 780nm.
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneAdjustLuminance.m' into MATLAB's
%      command window.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    01/02/18  jnm  Formatting & add checks for input variables
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
	scene = vcGetObject('scene');
	scene = sceneAdjustLuminance(scene, 100);  % Set to 100 cd/m2.
%}

% Check that required variables have been provided.
if notDefined('scene'), error('Error, scene is required.'); end
if notDefined('meanL'), error('Error, meanL is required.'); end

% Verify that current luminance exists, or calculate it
currentMeanL = sceneGet(scene, 'mean luminance');
photons = sceneGet(scene, 'photons');

photons = photons * (meanL / currentMeanL);

% We scale the photons and the illuminant data by the same amount to keep
% the reflectances in 0, 1 range.
scene = sceneSet(scene, 'photons', photons);
illuminant = sceneGet(scene, 'illuminant photons');
illuminant = illuminant * (meanL / currentMeanL);
scene = sceneSet(scene, 'illuminant photons', illuminant);

end