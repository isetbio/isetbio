% Demo off-axis  @cMosaic object computation
%
% Description:
%    Shows how to generate and use the new cone mosaic class, @cMosaic.
%    Here, we generate off-axis (located at different eccentricities)
%    cMosaic objects and compute their mean responses to a static stimulus
%    (rings/rays) with no eye movements.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxisDistortion
%   t_cMosaicBenchMark

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate the ring rays stimulus
fovDegs = 2;
scenePixels = 200;
wavefrontPixels = 201;
scene = sceneCreate('ringsrays', 10, scenePixels);
scene = sceneSet(scene, 'fov', fovDegs);
sceneWindow(scene);

%% Some optics parameters
turnOffLca = true;
diffractionLimitedHumanEye = false;
pupilSizeMm = 3;

%% Generate at a specified retinal location
%
% First coodinate is horizontal, second is vertical
mosaicEccDegs = [-2 0];

% Generate mosaic centered at target eccentricity
cm = cMosaic(...
    'sizeDegs', [1 1]*fovDegs, ...         % Mosaic size in degrees
    'eccentricityDegs', mosaicEccDegs ...  % Mosaic location in degrees
    );

%% Generate optics at same position
%
% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';   

% Optics database: choose between {'Polans2015', and 'Artal2012'}
% 
% The Polans2015 database has 10 subjects with measurements in a grid
% across the retina.
%
% The Artal2012 database has 41 subjects but only measured along the
% horizontal meridian.
opticsZernikeCoefficientsDataBase = 'Polans2015';            

% Select ranking of subject in database whose optics we will use.
% 1 is the best
subjectRankOrder = 1;

% Pick out the right subject from the database, and some info
% needed determined through our examination of the data.
switch (opticsZernikeCoefficientsDataBase)
    case 'Polans2015'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = PolansOptics.constants.subjectRanking;
        testSubjectID = rankedSujectIDs(subjectRankOrder);

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    case 'Artal2012'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
        testSubjectID = rankedSujectIDs(subjectRankOrder);

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
end

% Generate the optics for the right retinal location.
%
% It is the cMosaic object that knows how to process the databases
% and generate the optics.
[oiEnsemble, psfEnsemble] = ...
    cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
    'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
    'subjectID', testSubjectID, ...
    'pupilDiameterMM', pupilSizeMm, ...
    'zeroCenterPSF', false, ...
    'subtractCentralRefraction', subtractCentralRefraction, ...
    'wavefrontSpatialSamples', wavefrontPixels);
oi = oiEnsemble{1};
thePSFData = psfEnsemble{1};

% Get and manipulate the underlying wavefront data
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

% Adjust the pupil function to be diffraction limited.  You might
% want to do this if you were modeling an adaptive optics experiment.
if (diffractionLimitedHumanEye)
    zcoeffs = zeros(size(zcoeffs));
end
zcoeffs = wvfSet(wvf,'zcoeffs',zcoeffs);

% Put back the wvf into the oi. You need to compute on it first, though,
% if you actually changed anything about it.
wvf = wvfCompute(wvf);
wvfPlot(wvf,'psf','unit','um','wave',400,'plot range',40);
wvfPlot(wvf,'psf','unit','um','wave',550,'plot range',40);
wvfPlot(wvf,'psf','unit','um','wave',700,'plot range',40);
oi = oiSet(oi,'optics wvf',wvf);

% Compute the optical image of the scene
oi = oiCompute(oi,scene,'pad value','mean');
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



