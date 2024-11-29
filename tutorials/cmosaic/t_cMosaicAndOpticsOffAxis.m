% Demo off-axis  @cMosaic optics and cone excitation computation
%
% Description:
%    Shows how to generate and use the new cone mosaic class, @cMosaic.
%    Here, we generate off-axis (located at different eccentricities)
%    cMosaic objects and compute their mean cone excitations to a static stimulus
%    with no eye movements. In particular, this illustrates how to match up
%    and manipulate the optics as well.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxisDistortion
%   t_cMosaicBenchMark

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.
%    09/29/24  dhb  Simplified and extended for video.

%% Initialize
ieInit;
clear;
close all;

%% Generate a checkerboard with checks of specified size 
fovDegs = 1;
degreesPerCheck = 0.05;
nPixels = 200;
pixelsPerDegree = nPixels/fovDegs;
pixelsPerCheck = round(degreesPerCheck*pixelsPerDegree);
nCheckPairs = round(nPixels/(2*pixelsPerCheck));
scene = sceneCreate('checkerboard',pixelsPerCheck,nCheckPairs);
scene = sceneSet(scene, 'fov', fovDegs);
sceneWindow(scene);

%% Mosaic position on retina
%
% First coodinate is horizontal, second is vertical
mosaicEccDegs = [0 0];

%% Some optics parameters
turnOffLca = false;
diffractionLimitedHumanEye = false;
pupilDiamMM = 3;
addedDefocusDiopters = 0;

% Optics database: choose between {'Polans2015', 'Artal2012', 'Thibos2002'}
% 
% The Polans2015 database has 10 subjects with measurements in a grid
% across the retina.
%
% The Artal2012 database has 41 subjects but only measured along the
% horizontal meridian.
%
% Our curation of the Thibos2002 database has 70 subjects but only at the fovea.
opticsZernikeCoefficientsDataBase = 'Polans2015';  

% Rank order of subject's optics in database used.
% 1 is the best.  If you specify a number than the number of available
% subjects, the worst subject is used.
subjectRankOrder = 5;

%% Get optics information
%
% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';             

% Pick out the right subject from the database, and some info
% needed determined through our examination of the data.
switch (opticsZernikeCoefficientsDataBase)
    case 'Polans2015'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = PolansOptics.constants.subjectRanking;
        if (subjectRankOrder > length(rankedSujectIDs))
            testSubjectID = rankedSubjectIDs(end);
        else
            testSubjectID = rankedSujectIDs(subjectRankOrder);
        end

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    case 'Artal2012'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
        if (subjectRankOrder > length(rankedSujectIDs))
            testSubjectID = rankedSubjectIDs(end);
        else
            testSubjectID = rankedSujectIDs(subjectRankOrder);
        end

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
   
    case 'Thibos2002'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = ThibosOptics.constants.subjectRanking(whichEye);
        if (subjectRankOrder > length(rankedSujectIDs))
            testSubjectID = rankedSubjectIDs(end);
        else
            testSubjectID = rankedSujectIDs(subjectRankOrder);
        end

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ThibosOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
end

%% Generate mosaic centered at target eccentricity
cm = cMosaic(...
    'sizeDegs', [1 1]*fovDegs, ...         % Mosaic size in degrees
    'eccentricityDegs', mosaicEccDegs ...  % Mosaic location in degrees
    );

%% Generate the optics for the right retinal location.
%
% It is the cMosaic object that knows how to process the databases
% and generate the optics. The oi comes back in a cell array, and below
% we pick out the first and only element of that array. 
%
% The reason the oi comes back as a cell array is that it is possible to
% pass an eye movement path to this call, and then you need a shifted oi
% relative to the mosaic for every frame of the eye movement.  We will not
% explore that in this tutorial.
wavefrontPixels = 201;
oiEnsemble = ...
    cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
    'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
    'subjectID', testSubjectID, ...
    'pupilDiameterMM', pupilDiamMM, ...
    'zeroCenterPSF', true, ...
    'subtractCentralRefraction', subtractCentralRefraction, ...
    'wavefrontSpatialSamples', wavefrontPixels);
oi = oiEnsemble{1};

% Compute retinal magnification factor from focal length of the OI
% and check that it is matched to the cMosaic
focalLengthMeters = opticsGet(oiGet(oi,'optics'),'focal length');
micronsPerDegree = focalLengthMeters*tand(1)*1e6;
if (abs(micronsPerDegree-cm.micronsPerDegree)/micronsPerDegree > 1e-4)
    error('Mismatch between oi and cMosaic eye size')
end

%% Get and manipulate the underlying wavefront data
wvf = oiGet(oi,'optics wvf');
zcoeffs = wvfGet(wvf,'zcoeffs');

% Look at the current lca method, then set if desired.
% Options are 'humanlca', and 'none'.  You might want
% to set this to 'none' if you were modeling an experiment
% done in a way that compensated for human lca.
lcaMethod = wvfGet(wvf,'lcamethod');
if (turnOffLca)
    wvf = wvfSet(wvf,'lcamethod','none');
end

% Adjust the pupil function to be diffraction limited, if desired. You might
% want to do this if you were modeling an adaptive optics experiment.
if (diffractionLimitedHumanEye)
    zcoeffs = zeros(size(zcoeffs));
end

% Add in refractive error
addedDefocusMicrons = wvfDefocusDioptersToMicrons(addedDefocusDiopters, pupilDiamMM);
zcoeffs(5) = zcoeffs(5) + addedDefocusMicrons;
wvf = wvfSet(wvf,'zcoeffs',zcoeffs);

% Have a look at the point spread functions we are going to sue
wvf = wvfCompute(wvf);
wvfPlot(wvf,'psf','unit','um','wave',400,'plot range',40);
wvfPlot(wvf,'psf','unit','um','wave',550,'plot range',40);
wvfPlot(wvf,'psf','unit','um','wave',700,'plot range',40);

% Put the wvf back into the oi
oi = oiSet(oi,'optics wvf',wvf);

%% Compute the optical image of the scene
oi = oiCompute(oi,scene,'pad value','mean','crop',true);
oiWindow(oi);

%% Compute the noise-free excitation response
noiseFreeExcitationResponse = cm.compute(oi, 'opticalImagePositionDegs', mosaicEccDegs);

% Visualize mosaic response
cm.visualize('figureHandle', [], ...
    'axesHandle', [], ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'visualizedConeAperture', 'geometricArea', ...
    'domainVisualizationTicks', struct('x', -5:0.1:5, 'y', -5:0.1:5), ...
    'activation', noiseFreeExcitationResponse, ...
    'backgroundColor', 0.2*[1 1 1], ...
    'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEccDegs(1), mosaicEccDegs(2)));
