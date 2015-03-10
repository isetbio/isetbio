%% Script for testing the plotOI routine

%% 
s_initISET

%% Initialize the oi structure
scene = sceneCreate; 
scene = sceneSet(scene,'fov',4);
oi = oiCreate; oi = oiCompute(oi,scene);

%%
[uData, g] = plotOI(oi,'vline',[20 20]);

%%
[uData, g] = plotOI(oi,'hline',[20 20]);

%%
[uData, g] = plotOI(oi,'illuminance hline',[20 20]);

rows = round(oiGet(oi,'rows')/2);

%%
uData = plotOI(oi,' irradiance hline',[1,rows]);

%%
uData = plotOI(oi,'illuminance fft hline',[1,rows]);

%%
uData = plotOI(oi,'contrast hline',[1,rows]);

%%
uData = plotOI(oi,'irradiance image with grid',[],40);

%%
uData = plotOI(oi,'irradiance image wave',[],500,40);

%%
uData = plotOI(oi,'irradiance fft',[],450);

%%  Get some roiLocs
% uData = plotOI(oi,'irradiance energy roi');

%%
uData = plotOI(oi,'psf 550','um');

%%
uData = plotOI(oi,'otf 550','um');

%%
uData = plotOI(oi,'ls wavelength');

%% End