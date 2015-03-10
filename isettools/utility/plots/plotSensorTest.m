%% Script for testing the plotSensor routine

%% 
s_initISET

%% Initialize the sensor structure
scene = sceneCreate; 
scene = sceneSet(scene,'fov',4);
oi = oiCreate; oi = oiCompute(oi,scene);
vcAddAndSelectObject(oi);

sensor = sensorCreate; 
sensor = sensorSet(sensor,'qmethod', '10 bit');  % Linear, 10 bits
sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor);
sensorWindow('scale',1');

%%
plotSensor(sensor,'electrons hline',[20 20]);

%%
[uData, g] = plotSensor(sensor,'volts vline',[20 20]);

%%
[uData, g] = plotSensor(sensor,'volts hline',[20 20]);

%%
[uData, g] = plotSensor(sensor,'dv hline',[20 20]);

%%
[uData,g] = plotSensor(sensor,'sensor snr');

%% This requires hand intervention

uData = plotSensor(sensor,'electrons hist');
plotSensor(sensor,'electrons hist',uData.roiLocs)

%%
uData = plotSensor(sensor,'pixel snr');
%%
[uData, g] = plotSensor(sensor,'cfa block');
%%
[uData, g] = plotSensor(sensor,'cfa full');

%%
[uData, g] = plotSensor(sensor,'etendue');

%%
human = sensorCreate('human');
[uData, g] = plotSensor(human,'cone mosaic');

%% End


