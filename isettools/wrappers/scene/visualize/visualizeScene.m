function visualizeScene(scene, varargin)
p = inputParser;
p.addParameter('displayContrastProfiles', false, @islogical);
% Parse input
p.parse(varargin{:});
displayContrastProfiles = p.Results.displayContrastProfiles;


% retrieve the spatial support of the scene(in millimeters)
spatialSupportMilliMeters = sceneGet(scene, 'spatial support', 'mm');
% retrieve the sRGB components of the scene (just for visualization)
sceneRGBsettings = sceneGet(scene, 'rgb image');
% retrieve the XYZ tristimulus components of the scene
XYZmap = sceneGet(scene, 'xyz');
% compute the xy-chroma and the luminance maps
luminanceMap = squeeze(XYZmap(:,:,2));
xMap = squeeze(XYZmap(:,:,1))./sum(XYZmap,3);
yMap = squeeze(XYZmap(:,:,2))./sum(XYZmap,3);
% Compute mean luminance and mean chromaticity
meanLuminance = mean(luminanceMap(:));
meanChromaticity = [mean(xMap(:)) mean(yMap(:))];
% visualize the scene as RGB
visualizeSceneRGB(spatialSupportMilliMeters, 'mm', sceneRGBsettings, ...
    meanLuminance, meanChromaticity, sceneGet(scene, 'name'));

% retrieve the radiance of the scene as emitted photon rate 
% (photons/pixel/sec/nm)
scenePhotonRate = sceneGet(scene, 'photons');
% Retrieve wavelength support
wavelengthSupport = sceneGet(scene, 'wave');
wavelengthBandsToVisualize = 400:30:700;

% Visualize radiance maps
visualizeSceneRadiance(spatialSupportMilliMeters, 'mm', ...
    scenePhotonRate, wavelengthSupport, wavelengthBandsToVisualize);

if (displayContrastProfiles)
    visualizeSceneRadiance(spatialSupportMilliMeters, 'mm', ...
    scenePhotonRate, wavelengthSupport, wavelengthBandsToVisualize, ...
    'contrastProfilesOnly', true); 
end
