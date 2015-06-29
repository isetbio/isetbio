%% s_Toronto_Scenes

talkD = fullfile(isetbioRootPath,'scripts','toronto');

%% Scene section

% Show examples of test images
scene = sceneCreate('sweep frequency');
sceneSaveImage(scene,fullfile(talkD,'sweepFrequency.tif'),0.6);
scene = sceneCreate('freqorient');
sceneSaveImage(scene,fullfile(talkD,'images','frequencyOrient.tif'),0.6);

% Show examples of measured spectral images

% Show example of synthetic images created via Maya and RenderToolbox

% Illustrate  sceneFromFile with calibrated RGB data. 

%% Illustrate changing scene illuminant
sceneFile = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
scene = sceneFromFile(sceneFile,'multispectral');
scene = sceneAdjustLuminance(scene,50); % This sets the mean scene luminance
scene = sceneSet(scene,'fov',26.5); % match the scene field of view (fov) with the sensor fov

vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
scenePlot(scene,'illuminantEnergy')

%% Set illuminant to equal energy
scene = sceneAdjustIlluminant(scene,'equalEnergy.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow

%% Set to equal photon
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow

%% Convert the scene to the sunset color, Horizon_Gretag

% Notice that in this case 'Horizon_Gretag.mat' is a file name, not a
% data vector. 
scene = sceneAdjustIlluminant(scene,'Horizon_Gretag.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; 

