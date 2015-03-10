% s_sensorNoise
%
% Set sensor noise conditions and evaluate multiple samples.

%%
s_initISET

%%
scene = sceneCreate('slanted bar');
scene = sceneSet(scene,'fov',4);
vcAddAndSelectObject(scene); sceneWindow

oi = oiCreate;
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow

sensor = sensorCreate;

%% RE-WRITE using noise flag and reuse noise flag!@!!!
%
% This is out of date.
sensor = sensorSet(sensor,'exp time',0.050);
sensorNF = sensorComputeNoiseFree(sensor,oi);
v = sensorGet(sensorNF,'volts');

vcAddAndSelectObject(sensorNF); sensorImageWindow
% sensorNF = vcGetObject('sensor');

% Note that the largest amount of time is in the iePoisson calculator.
% Get a faster Poisson noise generator.
% profile on
nSamp = 100;
voltImages = sensorComputeSamples(sensorNF,nSamp);
% profile off
% profile viewer

%% Look at the noise histogram across all images
noiseImages = voltImages - repmat(v,[1 1 nSamp]);
vcNewGraphWin; hist(noiseImages(:),100)

%%
s = std(voltImages,0,3);
vcNewGraphWin;
imagesc(s); colorbar

%%
meanImage = mean(voltImages,3);
plot(meanImage(:),v(:),'.')
grid on;
axis equal

%%
s1 = sensorCompute(sensor,oi);
s2 = sensorCompute(sensor,oi);

v1 = sensorGet(s1,'volts');
v2 = sensorGet(s2,'volts');
hist(v1(:) - v2(:),100)




