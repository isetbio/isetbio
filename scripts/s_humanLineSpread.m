%% s_humanLineSpread
%
% Illustrate the photon absorptions from a line or grid stimulus. This
% is the practical measurement of the human line spread function,
% including photon noise, irregular spacing, and so forth.
%
% We illustrate the result for various spectral power distributions of the
% line stimuli, as well as the sum of a few nearby lines.
%
%

%% Initialize
ieInit

%% Create a line scene, human optics, and a human sensor
% This is a broad band stimulus, with a spectral power distribution of
% daylight, 6500 K.  We set the field of view to one degree.
%
% The optics are the estimated human optics, as per Marimont and Wandell in
% the mid-90s.
%
% The sensor has an approximation to the human cones, with random positions
% of the three cone types.

lineS = sceneCreate('line d65');
lineS = sceneSet(lineS,'h fov',1);

oi = oiCreate('wvf human');
oi = oiCompute(oi,lineS);

cm = coneMosaicRect;
cm.setSizeToFOV(1);
cm.integrationTime = 1;
cm.compute(oi);
cm.plot('hline absorptions lms','roi',[1 148]);

%% Change the line scene to 450nm and plot
preserveLuminance = 1;

line450S = sceneInterpolateW(lineS,450,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,line450S);
cm.compute(oi);
cm.plot('hline absorptions lms','roi',[1 148]);

%% Show the spread of the line at 550 nm
preserveLuminance = 1;

line550S = sceneInterpolateW(lineS,550,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,line550S);
cm.compute(oi);
cm.plot('hline absorptions lms','roi',[1 145]);
cm.plot('hline absorptions lms','roi',[1 146]);

%% Create a small grid pattern and image it on the sensor

imgSize = 128; lineSeparation = 32;
gridS = sceneCreate('gridlines',imgSize,lineSeparation);
gridS = sceneSet(gridS,'h fov',1);

oi = oiCreate('wvf human');
oi = oiCompute(oi,gridS);
oi = oiCrop(oi,'border');

cm.compute(oi);
cm.window;
cm.plot('hline absorptions lms','roi',[1 145]);

%% Show the 450 nm version of the grid.  Surprising, hunh?
% This version is very blurred, of course.  Surprisingly so.

grid450S = sceneInterpolateW(gridS,450,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,grid450S);
cm.compute(oi);
cm.window;
cm.plot('hline absorptions lms','roi',[1 145]);

%{
sensor = sensorSet(sensor,'exp time',0.050);
sensor = sensorCompute(sensor,oi);
sensor = sensorSet(sensor,'name','human-grid');

sensorWindow(sensor);
sensorPlot(sensor,'electrons hline',xy);
%}

%% END
