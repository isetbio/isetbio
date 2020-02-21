%%t_hyperspectralSceneTutorial  Illustrate ISETBio processing of hyperspectral scenes
%
% Description:
%     Illustrates how to:
%     - load a fetch a hyperspectral image
%     - generate a corresponding ISETbio scene
%     - generate a corresponding retinal optical image using wavefront-based
%       optics from a real subject.
%     - generate a realistic cone mosaic
%     - compute cone excitations for the given mosaic and scene to either
%       a statically viewed scene, or dynamically viewed scene (i.e. in the
%       presence of fixational eye movements.
%
%     This tutorial generates a number of figure which visualize different 
%     aspects of the processing pipeline, such as the spatial-spectral 
%     radiance of the scene and of the retinal image, the multi-spectral
%     point spread function of wavefront-based human optics, the cone mosaic 
%     spatial strucure, and the (dynamic) excitation of the cone mosaic.
%
%     
%  NPC, ISETBIO Team, 2016
%
% 07/24/18  npc  Wrote it.
% 10/08/18  dhb  Save cached mosaic in tempdir, not inside isetbio tree.


%% Close all open figures
% close all

%% Add subfolders to the path so we can have access to the supporting functions
% addpath(genpath(pwd));

%%
rdt = RdtClient('isetbio');
rdt.crp('/resources/scenes/hyperspectral/manchester_database/2004');
scene = rdt.readArtifact('scene3_2004_iset');
theScene = sceneFromBasis(scene);
%% Select scene #3 from the Manchester 2004 hyperspectral image dataset
sceneNo = 3;
sceneDataBase = 'manchester_database/2004';
loadSceneFromDisk = false;
if (loadSceneFromDisk)
    theScene = fetchScene(sceneNo, ...
        'origin', './resources', 'database', sceneDataBase);
else
    theScene = fetchScene(sceneNo, ...
        'origin', 'remote', 'database', sceneDataBase, ...
        'destination', '');
end

    
%% Parameters for the human optics
umPerDegree=300;

%% Generate human optics
theOI = oiCreate('wvf human');
    
%% Visualize the PSFs at some selected wavebands
visualizedWaveBands = [450 480 510 530 550 580 610 640 670];
visualizePSFfromOI(theOI, umPerDegree, ...
                'colormapToUse', gray(1024), ...
                'visualizedWavelengths', visualizedWaveBands, ...
                'rows', 3, 'cols', 3, ...
                'labelLastPSF', false, ...
                'displayWavelengthInTitle', ~false);
            
%% Compute the retinal image of the selected scene
theOI = oiCompute(theOI, theScene);

%% Visualize RGB renditions of the scene and of the retinal image
figNo = 1;
visualizedXRangeDegs = [-5 5]; visualizedYRangeDegs = [-3 3];
visualizeRGBrenditionsOfSceneAndRetinalImage(theScene, theOI, ...
    visualizedXRangeDegs, visualizedYRangeDegs, figNo);

%% Visualize the spatial spectral radiance of the scene at the selected wavelength bands
figNo = 3;
visualizeSpatialSpectralRadianceOfScene(theScene, visualizedWaveBands, ...
    visualizedXRangeDegs, visualizedYRangeDegs, figNo);

%% Visualize the spatial spectral radiance of the retinal image at the selected wavelength bands
figNo = 4;
visualizeSpatialSpectralRadianceOfRetinalImage(theOI, visualizedWaveBands, ...
    visualizedXRangeDegs, visualizedYRangeDegs, figNo);

%% Load a previously generated or generate a realistic cone mosaic
%% We should generate a database of high quality mosaics at some point
regenerateMosaic = true;
if (regenerateMosaic)
    mosaicFOV = 1.0;
    resamplingFactor = 7;  % 11 for high quality mosaic
    iterations = 10;  % This should be > 1000 for a high quality mosaic
    tol1 = 0.5;   % This should be < 0.1 for a high quality mosaic
    tol2 = 0.05;  % This should be < 0.01 for a high quality mosaic
    theConeMosaic = coneMosaicHex(resamplingFactor, ...
            'fovDegs', mosaicFOV*[1 0.5], ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', true, ...
            'latticeAdjustmentPositionalToleranceF', tol1, ...         
            'latticeAdjustmentDelaunayToleranceF', tol2, ... 
            'maxGridAdjustmentIterations', iterations);  
    % 5 millisecond integration time
    theConeMosaic.integrationTime = 5/1000;
    % Export mosaic
    save(fullfile(tempdir,'theMosaic'),  'theConeMosaic');
else
    load(fullfile(tempdir,'theMosaic'), 'theConeMosaic');
end

%% Visualize cone mosaic
theConeMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea');

%% Specify how many eye movements to generate. Specify 1 for no eye movements.
eyeMovementsNum = 50;
    
%% Generate eye movement path
if (eyeMovementsNum > 1)
    % Zero-centered fixational eye movement path
    fixEMobj = fixationalEM();
    fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsNum, ...
        'nTrials', 1, ...
        'rSeed', 857);
    emPath = fixEMobj.emPos;
else
    % Zero eye movements
    emPath = zeros(1,2,2);
end

%% Eye ofsset so we are centered on an interesting part of the image
eyeMovementOffsetDegs = [-1.5 0.4];
coneApertureMicrons = theConeMosaic.patternSampleSize(1)*1e6;
eyeMovementOffset = reshape(round(theConeMosaic.micronsPerDegree*eyeMovementOffsetDegs/coneApertureMicrons), [1 1 2]);
emPath = bsxfun(@plus, emPath, eyeMovementOffset);

%% Compute cone excitations for the current eye movement path
theConeExcitations = theConeMosaic.compute(theOI, 'emPath', emPath);

%% Visualize cone excitations response
figNo = 4;
visualizeConeExcitationResponse(theConeMosaic, theConeExcitations, ...
    theOI, emPath, figNo);

%%