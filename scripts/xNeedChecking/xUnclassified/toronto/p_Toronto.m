% Illustrate the spatial distribution of cone absorptions - somehow.
% This method interpolates the receptors (linearly) into the processor
% window.
% It should be possible to show them as an image in the sensor plane, too,
% I gues.  That would be new code.

%% 
s_initISET
talkD = fullfile(isetbioRootPath,'scripts','toronto');

% C:\Users\Wandell\Documents\Talks\20110715 Toronto\scripts
%%
p_Toronto_scene

p_Toronto_optics


%% Sensor section

% Illustrate cone mosaic and photon count

%% Illustrate chromatic aberration effects

scene = sceneCreate('gridlines');
scene = sceneSet(scene,'fov',1);
vcAddAndSelectObject(scene); sceneWindow;

oi = oiCreate('human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;

pSize = 2.5e-6;
sensor = sensorCreate('ideal',[],pSize,'human','bayer');
sensor = sensorSet(sensor,'exposure time',0.010);
sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor); sensorImageWindow;

%%  Cone absorption variance as a function of illuminant level or duration

% Illustrate the separation in the LMS clouds for pairs of wavelengths
oi = oiCreate('human');
pSize = 2.5e-6;
sensor = sensorCreate('ideal',[],pSize,'human','bayer');
sensor = sensorSet(sensor,'exposure time',0.010);

% This part is simple
scene = sceneCreate('uniform monochromatic',550);
oi = oiCompute(scene,oi);
sensor = sensorCompute(sensor,oi);

vcAddAndSelectObject(sensor); sensorImageWindow;

%
L = sensorGet(sensor,'electrons',1);
M = sensorGet(sensor,'electrons',2);
S = sensorGet(sensor,'electrons',3);
n = min(500,length(L));

% Make the same length
S = S(1:n); M = M(1:n); L = L(1:n);

vcNewGraphWin
plot3(L(:),M(:),S(:),'.')
xlabel('L-absorptions');
ylabel('M-Absorptions');
zlabel('S-absorptions');
axis square;
grid on



%% Implications for color matching

%% RGC ???  rgcVisualize ... with units

%% Statistics from surfaces to edges as a function of something.

% Light levels or what 
%

%% EJ's work on how cones drive these things


%% Aliasing in the cone mosaic ....

% Don't pass through the lens and show aliasing, like interferometry.

% Pass through the lens and protect from aliasing


%% Use L or M cone isolating gratings


%% Adapt within cone class option - adaptation in general


%% Checkerboard illusion

scene = sceneFromFile(fullfile(talkD,'images','Adelson','Illusion-Checkerboard.png'),'rgb',100,'lcdExample.mat');
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
vcAddAndSelectObject(scene); sceneWindow; % display sceneWindow

scene = sceneSet(scene,'fov',5);
oi = oiCreate('human');
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor,5);
sensor = sensorSet(sensor,'expTime',0.2);

oi = oiCompute(scene,oi);
sensor = sensorCompute(sensor,oi);

vcAddAndSelectObject(sensor); sensorImageWindow;





