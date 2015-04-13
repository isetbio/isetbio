function [srgb, XYZ] = sensorDemosaicCones(sensor, method, nFrames)
%% Demosaic and render human cone photon absorptions
%    srgb = sensorDemosaicCones(sensor, [method], [nFrames])
%
%  Inputs:
%    sensor  - sensor structure with photon absorption rates computed
%    method  - algorithm to be used for interpolation. Now support
%              'nearest', 'linear' (default), 'natural', 'cubic', 'v4'
%    nFrames - number of frames to be rendered, should be no larger than
%              number of frames in sensorGet(sensor, 'photons'). Default 1
%
%  Output:
%    srgb    - rendered srgb image
%    XYZ     - rendered XYZ image
%
%  Notes:
%    1) this function can only be used for demosaicing human cone mosaic.
%       For bayer patterns in cameras, use ISET camera modules instead
%    2) this function will detect dichromacy by checking the cone mosaic.
%       For dichromatic observers, the rendered image will be the
%       tranformed image for trichromats. See xyz2lms for more details
%       about dichromatic color transformation
%    3) For monochrome cone mosaic, this function will just scale it a gray
%       scale image
%
%  Examples:
%    fov = 1;
%    scene = sceneCreate; scene = sceneSet(scene, 'h fov', fov);
%    oi = oiCreate('human'); oi = oiCompute(oi, scene);
%    sensor = sensorCreate('human'); sensor = sensorCompute(sensor, oi);
%    srgb = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);
%
% (HJ) ISETBIO TEAM, 2015