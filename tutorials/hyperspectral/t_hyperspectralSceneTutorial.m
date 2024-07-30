% t_hyperspectralSceneTutorial  
%
% Description:
%     Illustrates how to load a hyperspectral image, compute its optical
%     image and a cMosaic excitation response.
%     
%  NPC, ISETBIO Team, 2016
%
% 07/24/18  npc  Wrote it.
% 10/08/18  dhb  Save cached mosaic in tempdir, not inside isetbio tree.
% 01/02/24  npc  cMosaic version
% 05/18/24  bw   Reduced size for speed, removed dropbox directory.

%% Close all open figures
clear; close all

%% Load hyperspectral image data
defaultImage = false;
if (defaultImage)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
else
    % This will work if you are in the Brainard Lab and have the
    % HyperspectralSceneTutorial folder on your lab dropbox path
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, 'scene3.mat'), 'scene');
end

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

ieNewGraphWin;
image(xSupportDegs, ySupportDegs, sceneRGB);

% Generata a 3x3 cMosaic (which covers part of the image only)
theConeMosaic = cMosaic('sizeDegs', [1,1]);

% Compute its activation
theConeMosaicActivation = theConeMosaic.compute(theOI);

% Visualize the cone mosaic
theConeMosaic.visualize();

% Visualize its activation to the hyperspecral image
theConeMosaic.visualize(...
    'activation', theConeMosaicActivation);

%% END