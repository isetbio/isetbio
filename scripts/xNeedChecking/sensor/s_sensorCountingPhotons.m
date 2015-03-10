%% s_sensorCountingPhotons
%
% Let's explore how many scene photons find their way down to a pixel in a
% sensor.  Used for some teaching.
%
% (c) Imageval, 2011

%%
s_initISET

%%
sFile = fullfile(isetRootPath,'data','images','rgb','hats.jpg');
scene = sceneFromFile(sFile,'rgb');
scene = sceneAdjustIlluminant(scene,'D65.mat');
vcAddAndSelectObject(scene);
sceneWindow;

%% Plot a region of interest from the scene.  

% You can do the plot and then get the roiRect from the get(gcf,'userdata')

% This region is in the middle of the body of the animal
roiRect = [225 164 16 16];
roiLocs = ieRoi2Locs(roiRect);
[udata, f] = plotScene(scene,'radiance photons roi',roiLocs);

% The sum of the mean number of photons from all the wavelengths 
% q/s/sr/nm/m2
totalQ = sum(udata.photons(:))

%% Create the optical image

oi = oiCreate;
optics = oiGet(oi,'optics');
roiRect = [291 202 16 23];

% F number 2
optics = opticsSet(optics,'f number',2);
oi = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;

plotOI(oi,'irradiance photons roi',roiRect);
udata = get(gcf,'userdata');
totalQ = sum(udata.y);
fprintf('F# %f - Total Q: %e\n',opticsGet(optics,'f number'),totalQ)

% F number 4
optics = opticsSet(optics,'f number',4);
oi = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;
plotOI(oi,'irradiance photons roi',roiRect);
udata = get(gcf,'userdata');
totalQ = sum(udata.y);
fprintf('F# %f - Total Q: %e\n',opticsGet(optics,'f number'),totalQ)

% F number 8
optics = opticsSet(optics,'f number',8);
oi = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;
plotOI(oi,'irradiance photons roi',roiRect);
udata = get(gcf,'userdata');
totalQ = sum(udata.y);
fprintf('F# %f - Total Q: %e\n',opticsGet(optics,'f number'),totalQ)



