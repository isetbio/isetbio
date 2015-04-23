function t_GenerateConeMosaicImage()
%
% Show how to make a nice viewable image of a cone mosaic.
%
% 4/21/15 npc   Wrote it.
% 4/22/15 dhb   Tweaks for isetbio compatibility.
% 4/23/15 npc   Minor fixes.

%% Clear
close all; clear global; ieInit;

%% Create human sensor to get a cfa
params.rgbDensities = [0.0 0.625 0.325 .05];
sensor = sensorCreate('human',[],params);
sensor = sensorSet(sensor, 'noise flag', 0);
sensor = sensorSet(sensor,'exp time',2);
sensor = sensorSet(sensor,'rows',128);
sensor = sensorSet(sensor,'cols',128);
pixel = sensorGet(sensor,'pixel');
pixel = pixelSet(pixel,'sizesamefillfactor',(200/120)*[1.5e-6 1.5e-6]);
sensor = sensorSet(sensor,'pixel',pixel);
coneCFAPattern = sensorGet(sensor,'cfa pattern');

%% Extract a reasonable size region
cfaSize = 30;
coneCFAPattern = coneCFAPattern(1:cfaSize,1:cfaSize);

%% Specify cone aperture size
coneSize = 15;

%% generate the cone mosaic image
[coneMosaicStandardImage,coneCFARawImage] = generateConeMosaicImage(coneCFAPattern, coneSize, 'standard');
coneMosaicWRImage = generateConeMosaicImage(coneCFAPattern, coneSize, 'williams_roorda');

% Show the results
h = figure(1); set(h, 'Name', 'raw cfa mosaic', 'Position', [10 10 100 100]); clf
imshow(coneCFARawImage);

h = figure(2); set(h, 'Name', 'standard style cone mosaic', 'Position', [400 400 100 100]); clf
imshow(coneMosaicStandardImage); truesize;

h = figure(3); set(h, 'Name', 'Williams/Roorda style cone mosaic', 'Position', [860 400 100 100]); clf
imshow(coneMosaicWRImage); truesize;

end

