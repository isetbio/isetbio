function visualizeOpticalImage(opticalImage, varargin)

p = inputParser;
p.addParameter('displayRetinalContrastProfiles', false, @islogical);
p.addParameter('displayRadianceMaps', true, @islogical);
p.addParameter('axesHandle', []);
p.addParameter('crossHairsAtOrigin', false, @islogical);
p.addParameter('avoidAutomaticRGBscaling', false, @islogical);

% Parse input
p.parse(varargin{:});
displayRetinalContrastProfiles = p.Results.displayRetinalContrastProfiles;
displayRadianceMaps = p.Results.displayRadianceMaps;
crossHairsAtOrigin = p.Results.crossHairsAtOrigin;

% retrieve the spatial support of the scene(in millimeters)
spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
% Convert spatial support in degrees
optics = oiGet(opticalImage, 'optics');
focalLength = opticsGet(optics, 'focal length');
mmPerDegree = focalLength*tand(1)*1e3;
spatialSupportDegs = spatialSupportMM/mmPerDegree;


% retrieve the sRGB components of the optical image (just for visualization)
rgbImage = oiGet(opticalImage, 'rgb image');
% retrieve the XYZ tristimulus components of the optical image
XYZmap = oiGet(opticalImage, 'xyz');
xMap = squeeze(XYZmap(:,:,1))./sum(XYZmap,3);
yMap = squeeze(XYZmap(:,:,2))./sum(XYZmap,3);
% Compute mean chromaticity
meanChromaticity = [mean(xMap(:)) mean(yMap(:))];
optics = oiGet(opticalImage, 'optics');
opticsName = opticsGet(optics, 'name');
ax = visualizeSceneRGB(spatialSupportDegs, 'degs', rgbImage, ...
    [], meanChromaticity, ' ', 'axesHandle', p.Results.axesHandle, ...
    'avoidAutomaticRGBscaling', p.Results.avoidAutomaticRGBscaling);

if (crossHairsAtOrigin)
    % Obtain spatial support in mm
    spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
    % Compute mm/deg conversion factor
    optics = oiGet(opticalImage, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
      
    % Convert to degs
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportX = spatialSupportDegs(1,:,1);
    spatialSupportY = spatialSupportDegs(:,1,2);
    hold(ax, 'on');
    plot(ax,[spatialSupportX(1) spatialSupportX(end)], [0 0], 'k-');
    plot(ax,[0 0], [spatialSupportY(1) spatialSupportY(end)],  'k-');
end

if (displayRadianceMaps)
    % retrieve the retinal irradiance (photon rate in photons/pixel/sec/nm)
    retinalImagePhotonRate = oiGet(opticalImage, 'photons');

    % retrieve wavelength support
    wavelengthSupport = oiGet(opticalImage, 'wave');
    wavelengthBandsToVisualize = 400:30:700;

    % Visualize radiance maps
    visualizeSceneRadiance(spatialSupportDegs, 'degs', ...
        retinalImagePhotonRate, wavelengthSupport, wavelengthBandsToVisualize);
end

% Visualize radiance contrast profiles
if (displayRetinalContrastProfiles) && (displayRadianceMaps)
    visualizeSceneRadiance(spatialSupportDegs, 'degs', ...
        retinalImagePhotonRate, wavelengthSupport, wavelengthBandsToVisualize, ...
        'contrastProfilesOnly', true);     
end

end
