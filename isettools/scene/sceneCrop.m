function [scene, rect] = sceneCrop(scene, rect)
% Crop scene data.
%
% Syntax:
%	[scene, rect] = sceneCrop(scene, [rect])
%
% Description:
%    Crop the data (photons) in the scene or optical image to within the
%    rectangle, rect. If rect is not defined, then a graphical routine is
%    initiated for selecting the rectangle. Following that, the values of
%    rect can be returned.
%
%    Because these are multispectral data, we can't use the usual imcrop.
%    Instead, we use vcROISelect to return the selected data. Then we turn
%    the selected data into a photon image and reset the relevant
%    parameters in the scene structure.
%
%    N.B. The source contains executable examples of usage, which can be
%    accessed by typing 'edit sceneCrop.m' in the command window.
%
%    There are examples in the code. Type 'edit sceneCrop' into the Command
%    Window to access.
%
% Inputs:
%    scene - The scene structure
%    rect  - (Optional) The rectangle to crop the data/image within
%
% Outputs:
%    scene - The modified scene structure
%    rect  - The rectangle's information
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/22/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match the Wiki.

% Examples:
%{
    % ETTBSkip.  Requires user input, so do not auto run this example
    [val, scene] = vcGetSelectedObject('SCENE');
	newScene = sceneCrop(scene);
%}

if notDefined('scene'), error('You must define a scene.'); end

if notDefined('rect')
    [roiLocs, rect] = vcROISelect(scene);
else
    cmin = rect(1);
    cmax = rect(1) + rect(3);
    rmin = rect(2);
    rmax = rect(2) + rect(4);
    [c, r] = meshgrid(cmin:cmax, rmin:rmax);
    roiLocs = [r(:), c(:)];
end

% The number of selected columns and rows
c = rect(3) + 1;
r = rect(4) + 1;
% wave = sceneGet(scene, 'nwave');

% These are in XW format.
photons = vcGetROIData(scene, roiLocs, 'photons');
photons = XW2RGBFormat(photons, r, c);

% Now build up the new object.
scene = sceneClearData(scene);
scene = sceneSet(scene, 'photons', photons);
[luminance, meanL] = sceneCalculateLuminance(scene);
scene = sceneSet(scene, 'luminance', luminance);
scene = sceneSet(scene, 'meanLuminance', meanL);

end