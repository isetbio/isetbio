% s_opticsDefocusScene
%
% Illustrate how to apply defocus to a scene.  
%
% The first section shows how to create a defocused image when the defocus
% is specified in diopters.
%
% The next section shows how to calculate the blur caused by an image that
% is mis-focused from the sensor plane by a specific distance.
%
% Copyright ImagEval Consultants, LLC, 2013

%%
s_initISET
wbStatus = ieSessionGet('waitbar');
ieSessionSet('waitbar','on');

%% Create a test scene

% This is a simple picture that sweeps out 5 deg of visual angle
wave = 400:10:700;
fullFileName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs');
scene = sceneFromFile(fullFileName,'multispectral',[],[],wave);
scene = sceneSet(scene,'fov',5);

maxSF = sceneGet(scene,'max freq res','cpd');
nSteps = min(ceil(maxSF),70);              % Round up, but don't go too high.
sampleSF = linspace(0, maxSF, nSteps);     % cyc/deg
vcAddAndSelectObject(scene); sceneWindow;

%% Standard shift invariant optics
% We are assuming a diffraction limited lens with some defocus.
%
oi     = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'model','shift invariant');
wave   = opticsGet(optics,'wave');

%% Create the optical transfer function (OTF) for the specified defocus 

% Initialize the defocus for each wavelength
defocus = zeros(size(wave));
D = 5;                     % Defocus
defocus = defocus + D;     % In units of diopters

% Create the defocused otf 
[otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,defocus);

% This is the OTF as a function of wavelength
vcNewGraphWin; mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%% Calculate the optical image given the defocus

% First create the optics with this defocus OTF
optics = opticsBuild2Dotf(optics,otf,sampleSFmm);

oi = oiSet(oi,'optics',optics);
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name',sprintf('%.1f-defocus',D));
vcAddAndSelectObject(oi); oiWindow;

%% Perfect focus, and then two cases where the distance is wrong

% Here is the perfect image with no defocus.
defocus = zeros(size(wave));
D = 0;
defocus = defocus + D;
[otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,defocus);

optics = opticsBuild2Dotf(optics,otf,sampleSFmm);

oi = oiSet(oi,'optics',optics);
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name',sprintf('%.1f-defocus',D));

v1 = vcAddAndSelectObject(oi); oiWindow;
%% Calculate the consequences for an image that is 10 um off

% Correct focal length
fLength = opticsGet(optics,'focal length');
lensPower = opticsGet(optics,'diopters');

% Suppose the image focal length misses by this much 10 microns
deltaDistance = 10e-6;

% The effective power of a lens imaging at that distance is
actualPower = 1 / (fLength - deltaDistance);

% Correction needed - which is also the amount of defocus
D = actualPower - lensPower;
defocus = zeros(size(wave));
defocus = defocus + D;     % In units of diopters

% Create the defocused otf 
[otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,defocus);
optics = opticsBuild2Dotf(optics,otf,sampleSFmm);

% Compute and name and view
oi = oiSet(oi,'optics',optics);
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name',sprintf('%.1f-defocus',D));

v2 = vcAddAndSelectObject(oi); oiWindow;

%% Again, but for for an image that is 40 um off

% Correct focal length
fLength = opticsGet(optics,'focal length');
lensPower = opticsGet(optics,'diopters');

% Suppose the image focal length misses by this much 10 microns
deltaDistance = 40e-6;

% The effective power of a lens imaging at that distance is
actualPower = 1 / (fLength - deltaDistance);

% Correction needed - which is also the amount of defocus
D = actualPower - lensPower;
defocus = zeros(size(wave));
defocus = defocus + D;     % In units of diopters

% Create the defocused otf 
[otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,defocus);
optics = opticsBuild2Dotf(optics,otf,sampleSFmm);

% Compute and name and view
oi = oiSet(oi,'optics',optics);
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name',sprintf('%.1f-defocus',D));

v3 = vcAddAndSelectObject(oi); oiWindow;


%% Show the three optical image results in three windows
imageMultiview('oi',[v1,v2,v3],true);
ieSessionSet('waitbar',wbStatus);

%% End