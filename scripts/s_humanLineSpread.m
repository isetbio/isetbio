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
oi = oiCompute(oi,lineS,'pad value','mean');

cm = coneMosaicRect;
cm.setSizeToFOV(1);
cm.integrationTime = 1;
cm.compute(oi);
[~,target] = cm.plot('hline absorptions lms','roi',[1 148]);
title(target.Children(end),'Wave: broadband');


%% Change the line scene to 450nm and plot
preserveLuminance = 1;

line450S = sceneInterpolateW(lineS,450,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,line450S,'pad value','mean');
cm.compute(oi);
[~,target] = cm.plot('hline absorptions lms','roi',[1 148]);
title(target.Children(end),'Wave: 450nm');

%% Show the spread of the line at 550 nm
preserveLuminance = 1;

line550S = sceneInterpolateW(lineS,550,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,line550S,'pad value','mean');
cm.compute(oi);
cm.plot('hline absorptions lms','roi',[1 145]);
title('Wave: 550nm');

[~,target] = cm.plot('hline absorptions lms','roi',[1 146]);
title(target.Children(end),'Wave: 550nm');

%% Create a small grid pattern and image it on the sensor

imgSize = 128; lineSeparation = 32;
gridS = sceneCreate('gridlines',imgSize,lineSeparation);
gridS = sceneSet(gridS,'h fov',1);

oi = oiCreate('wvf human');
oi = oiCompute(oi,gridS,'pad value','mean');
oi = oiCrop(oi,'border');

cm.compute(oi);
cm.window;
[~,target] = cm.plot('hline absorptions lms','roi',[1 145]);
title(target.Children(end),'Wave: broadband');

%% Show the 450 nm version of the grid.  Surprising, hunh?
% This version is very blurred, of course.  Surprisingly so.

grid450S = sceneInterpolateW(gridS,450,preserveLuminance);
oi = oiCreate('wvf human');
oi = oiCompute(oi,grid450S,'pad value','mean');
cm.compute(oi);
cm.window;
[~,target] = cm.plot('hline absorptions lms','roi',[1 145]);
title(target.Children(end),'Wave: 450nm');

%% END
