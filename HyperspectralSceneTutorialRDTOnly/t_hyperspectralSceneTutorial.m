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
        'destination', './resources');
end

    theScene = fetchScene(sceneNo, ...
        'origin', 'remote', 'database', sceneDataBase, ...
        'destination', './resources');
    
%% Parameters for the human optics
oiParams = struct(...
    'opticsModel', 'ThibosDefaultSubject3MMPupil', ...
    'wavefrontSpatialSamples', 261*2+1, ...
    'pupilDiamMm', 3.0, ...
    'umPerDegree', 300);

%% Generate human optics
theOI = oiWithCustomOptics(oiParams.opticsModel, ...
    oiParams.wavefrontSpatialSamples, ...
    oiParams.pupilDiamMm, ...
    oiParams.umPerDegree);
    
%% Visualize the PSFs at some selected wavebands
visualizedWaveBands = [450 480 510 530 550 580 610 640 670];
visualizePSFfromOI(theOI, oiParams.umPerDegree, ...
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
regenerateMosaic = ~true;
if (regenerateMosaic)
    mosaicFOV = 2.0;
    theConeMosaic = coneMosaicHex(11, ...
            'fovDegs', mosaicFOV*[1 0.5], ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', true, ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...         
            'latticeAdjustmentDelaunayToleranceF', 0.05, ... 
            'maxGridAdjustmentIterations', 500);
    % 5 millisecond integration time
    theConeMosaic.integrationTime = 5/1000;
    % Export mosaic
    save('theMosaic',  'theConeMosaic');
else
    load('theMosaic', 'theConeMosaic');
end

%% Visualize cone mosaic
theConeMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea');

%% Specify how many eye movements to generate. Specify 1 for no eye movements.
eyeMovementsNum = 20;
    
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