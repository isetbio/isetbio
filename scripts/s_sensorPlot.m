%% Script for testing the sensorPlot routine

%% 
ieInit;

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
sensorPlot(sensor,'photons hline',[20 20]);

%%
[uData, g] = sensorPlot(sensor,'volts vline',[20 20]);

%%
[uData, g] = sensorPlot(sensor,'volts hline',[20 20]);

%%
[uData, g] = sensorPlot(sensor,'dv hline',[20 20]);

%%
[uData,g] = sensorPlot(sensor,'sensor snr');

%% This requires hand intervention

uData = sensorPlot(sensor,'photons hist');
sensorPlot(sensor,'photons hist',uData.roiLocs)

%%
uData = sensorPlot(sensor,'pixel snr');
%%
[uData, g] = sensorPlot(sensor,'cfa block');
%%
[uData, g] = sensorPlot(sensor,'cfa full');

%%
[uData, g] = sensorPlot(sensor,'etendue');

%%
human = sensorCreate('human');
[uData, g] = sensorPlot(human,'cone mosaic');

%% End


