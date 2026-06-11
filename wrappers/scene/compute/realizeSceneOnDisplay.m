function realizedScene = realizeSceneOnDisplay(scene, display)
% Method to generate an ISETBio scene as realized on the passed display
%
% Syntax:
%   realizedScene = realizeSceneOnDisplay(scene, display)
%
% Description:
%    This function uses an idealized ISETBio scene, in which radiance can be 
%    specified independently at each wavelength, and generates a metameric
%    ISETBio scene that can be realized on the passed display. The two
%    scenes have the same luminance and chromaticity, but differ in the
%    spectral radiance composition, as the latter scene is generated based
%    on the spectral power distributions of the passed display's primaries.
%
% Inputs:
%   scene:        - The input ISETBio scene 
%   display:      - The display on which the input scene is to be realized
%
% Outputs:
%   realizedScene - The ISETBio scene corresponding to the realized scene
%                   on the passed display.
%
% History:
%    11/23/18  NPC  ISETBIO TEAM, 2018

% Extract the input scene's XYZ components
sceneXYZ = sceneGet(scene, 'xyz');
    
% Extract the display's RGB-to-XYZ transformation matrix
displayRGBtoXYZ = displayGet(display, 'rgb2xyz');
    
% Generate linear RGB primaries necessary to reproduce the scene's XYZ component.
sceneRGBPrimaries = imageLinearTransform(sceneXYZ, inv(displayRGBtoXYZ));

% Check for any rgbSettings values above 1.0 and issue a warning that some
%values will be clipped to 1.0
if (any(sceneRGBPrimaries(:)>1.0))
    warndlg('Image is out of gamut (> 1))','Clipping to gamut');
    sceneRGBPrimaries(sceneRGBPrimaries>1.0) = 1;
    fprintf('Primaries range: %2.2f - %2.2f\n', ...
    min(sceneRGBPrimaries(:)), max(sceneRGBPrimaries(:)));

end
    
if (any(sceneRGBPrimaries(:)<0.0))
    warndlg('Image is out of gamut (< 0))', 'Clipping to gamut');
    sceneRGBPrimaries(sceneRGBPrimaries<0.0) = 0;
    fprintf('Primaries range: %2.2f - %2.2f\n', ...
    min(sceneRGBPrimaries(:)), max(sceneRGBPrimaries(:)));
end

% Extract inverse gamma table 
inverseGammaTable = displayGet(display, 'inverse gamma');
% Normalize it
inverseGammaTable = inverseGammaTable/max(inverseGammaTable(:));

% Pass linear RGB primaries via inverse gamma to generate the display
% settings values
sceneSettings = ieLUTLinear(sceneRGBPrimaries, inverseGammaTable);
 
% Generate a scene based on these RGB settings
meanLuminance = [];
realizedScene = sceneFromFile(sceneSettings, 'rgb', meanLuminance, display);

% Set the realized scene size and view distance to match those of the scene
realizedScene = sceneSet(realizedScene, 'wangular', sceneGet(scene,'wangular'));
realizedScene = sceneSet(realizedScene, 'distance', sceneGet(scene,'distance'));
end
