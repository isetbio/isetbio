function visualizeSpatialSpectralRadianceOfScene(theScene, visualizedWaveBands, ...
            visualizedXRangeDegs, visualizedYRangeDegs, figNo)
% Visualize the spatial spectral radiance of a scene
%
% 7/25/18  npc  Wrote it
%
    spatialSizeDegs = [sceneGet(theScene, 'wAngular') sceneGet(theScene, 'hAngular') ]; 
    radianceWaveSupport = sceneGet(theScene, 'wave');
    radiancePhotons = sceneGet(theScene, 'photons');
    visualizeSpatialSpectralRadiance(spatialSizeDegs, radianceWaveSupport, radiancePhotons, visualizedWaveBands, visualizedXRangeDegs, visualizedYRangeDegs, figNo, 'scene');
end

