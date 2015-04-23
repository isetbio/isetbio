function t_GenerateConeMosaicImage()
%
% Show how to make a nice viewable image of a cone mosaic.
%
% 4/21/15 ncp   Wrote it.
% 4/22/15 dhb   Tweaks for isetbio compatibility.

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
coneSize = 16;

%% generate the cone mosaic image
[coneMosaicStandardImage,coneCFARawImage] = generateConeMosaicImage(coneCFAPattern, coneSize, 'standard');
coneMosaicWRImage = generateConeMosaicImage(coneCFAPattern, coneSize, 'williams_roorda');

% Show the results
h = figure(1); set(h, 'Name', 'raw cfa mosaic'); clf
imshow(coneCFARawImage);

h = figure(2); set(h, 'Name', 'standard style cone mosaic'); clf
imshow(coneMosaicStandardImage); truesize;

h = figure(3); set(h, 'Name', 'willaims/roorday style cone mosaic'); clf
imshow(coneMosaicWRImage); truesize;

end

