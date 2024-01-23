%%t_hyperspectralSceneTutorial  
%
% Description:
%     Illustrates how to load a hyperspectral image, compute its optical
%     image and a cMosaic excitation response.
%     
%  NPC, ISETBIO Team, 2016
%
% 07/24/18  npc  Wrote it.
% 10/08/18  dhb  Save cached mosaic in tempdir, not inside isetbio tree.
% 01/02/25  npc  cMosaic version

%% Close all open figures
% close all

%% Load hyperspectral image data from Dropbox
dropboxDirPath = localDropboxDir();
scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
load(fullfile(scenesDir, 'scene3.mat'), 'scene');

%% Generate human optics
theOI = oiCreate('wvf human');
              
%% Compute the retinal image of the selected scene
theOI = oiCompute(theOI,scene,'pad value','mean');

%% Visualize the scene
sceneRGB = sceneGet(scene, 'rgb');
dxy = sceneGet(scene, 'angular resolution');
[rows, cols, ~] = size(sceneRGB);
xSupportDegs = (1:cols)*dxy(1);
ySupportDegs = (1:rows)*dxy(2);
xSupportDegs = xSupportDegs - mean(xSupportDegs);
ySupportDegs = ySupportDegs - mean(ySupportDegs);
figure(); image(xSupportDegs, ySupportDegs, sceneRGB);

% Generata a 3x3 cMosaic (which covers part of the image only)
theConeMosaic = cMosaic('sizeDegs', [3,3]);

% Compute its activation
theConeMosaicActivation = theConeMosaic.compute(theOI);

% Visualize the cone mosaic
theConeMosaic.visualize();

% Visualize its activation to the hyperspecral image
theConeMosaic.visualize(...
    'activation', theConeMosaicActivation);