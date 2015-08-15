%% s_Toronto_Scenes

%%
talkD = 'C:\Users\Brian\Talks\20110715 Toronto';

%% Make tif images of test scenes

% Show examples of test images
scene = sceneCreate('sweep frequency');
sceneSaveImage(scene,fullfile(talkD,'images','sweepFrequency.tif'),0.6);

scene = sceneCreate('freqorient');
sceneSaveImage(scene,fullfile(talkD,'images','frequencyOrient.tif'),0.6);

scene = sceneCreate('gridLines');
sceneSaveImage(scene,fullfile(talkD,'images','gridLines.tif'),0.6);

scene = sceneCreate('slanted bar');
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
% vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
sceneSaveImage(scene,fullfile(talkD,'images','slantedBar.tif'),0.6);

scene = sceneCreate('checkerboard');
sceneSaveImage(scene,fullfile(talkD,'images','checkerboard.tif'),0.6);

scene = sceneCreate('macbethd65');
sceneSaveImage(scene,fullfile(talkD,'images','macbethd65.tif'),0.6);

scene = sceneCreate('point array');
sceneSaveImage(scene,fullfile(talkD,'images','pointArray.tif'),0.6);



%% Illustrate changing scene illuminant
sceneFile = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
scene = sceneFromFile(sceneFile,'multispectral');
scene = sceneAdjustLuminance(scene,50); % This sets the mean scene luminance
scene = sceneSet(scene,'fov',26.5); % match the scene field of view (fov) with the sensor fov

vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
scenePlot(scene,'illuminantEnergy')
% Show examples of measured spectral images
% ?? scenePlot(scene,'illuminantEnergy')

% Set illuminant to equal energy
scene = sceneAdjustIlluminant(scene,'equalEnergy.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
scenePlot(scene,'illuminantEnergy')

% Set to equal photon
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
scenePlot(scene,'illuminantEnergy')

% Convert the scene to the sunset color, Horizon_Gretag
% Notice that in this case 'Horizon_Gretag.mat' is a file name, not a
% data vector. 
scene = sceneAdjustIlluminant(scene,'Horizon_Gretag.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; 
scenePlot(scene,'illuminantEnergy')

%% Show example of synthetic images created via Maya and RenderToolbox
fName = fullfile(isetRootPath,'data','images','multispectral','piano3d.mat');
load(fName);
scene = sceneSet(scene,'fov',3);
vcAddAndSelectObject('scene',scene); sceneWindow; 

scenePlot(scene,'depth map contour');

%% Illustrate sceneFromFile with  RGB data on a calibrated display
fName = fullfile(isetRootPath,'data','images','rgb','eagle.jpg');
scene = sceneFromFile(fName,'rgb',100,'lcdExample.mat');
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; 

sz = sceneGet(scene,'size');
xy = [1,round(sz(1)/2)];

scenePlot(scene,'hline radiance',xy);
scenePlot(scene,'hline luminance',xy);

scene = sceneSet(scene,'distance',100);
gridSpacing = 2000; %mm
scenePlot(scene,'radiance image with grid',[],gridSpacing);

