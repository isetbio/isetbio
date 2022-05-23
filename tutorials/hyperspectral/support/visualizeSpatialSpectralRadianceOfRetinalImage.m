function visualizeSpatialSpectralRadianceOfRetinalImage(theOI, visualizedWaveBands, ...
            visualizedXRangeDegs, visualizedYRangeDegs, figNo)
% Visualize the spatial spectral radiance of a retinal image
%
% 7/25/18  npc  Wrote it
%
    spatialSizeDegs = [oiGet(theOI, 'wAngular') oiGet(theOI, 'hAngular') ]; 
    radianceWaveSupport = oiGet(theOI, 'wave');
    radiancePhotons = oiGet(theOI, 'photons');
    visualizeSpatialSpectralRadiance(spatialSizeDegs, radianceWaveSupport, radiancePhotons, visualizedWaveBands, visualizedXRangeDegs, visualizedYRangeDegs, figNo, 'retinal image');
end