function fName = sceneSaveImage(scene, fName, gam)
% Write png image approximating appearance of photon data
%
% Syntax:
%	fullpathName = sceneSaveImage([scene], [fullpathname], [gam]);
%
% Description:
%    Save out an RGB image of the photon image as a png file. If the name
%    is not passed in, then the user is queried to select the fullpath name
%    of the output file. This routine is used for scenes. oiSaveImage is
%    used for optical images.
%
%    The rgb image is obtained from the scene via a call to sceneGet, which
%    gets the photons, converts them to an XYZ image, and then converts the
%    XYZ image to an sRGB format. This is the format that is displayed in
%    the scene window, and it is the sRGB format that is saved here.
%
% Inputs:
%    scene - (Optional) The scene structure. Default is retrieve a scene.
%    fName - (Optional) The filename. Default is to query user to select.
%    gam   - (Optional) The gamma. Default is to not provide.
%
% Outputs:
%    fName - The filename containing the saved scene
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/21/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match the Wiki.


if notDefined('scene'), scene = vcGetObject('scene'); end

% Get RGB file name (tif)
if notDefined('fName')
    fName = vcSelectDataFile('stayput', 'w', 'png', 'Image file (png)');
end

% open window.
if notDefined('gam')
    RGB = sceneGet(scene, 'rgb image');
else
    RGB = sceneGet(scene, 'rgb image', gam);
end

imwrite(RGB, fName, 'png');

end