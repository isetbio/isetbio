% Script to examine the properties of compute-ready mosaics 
% This script is also used to generate materials for the validation
% figures of our PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicCheckerBoard  
%}

% Initialize session
close all; clear all;

% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
%ASPPManager = AppleSiliconParPoolManager(12);
%
ASPPManager = AppleSiliconParPoolManager('half max');
%ASPPManager = AppleSiliconParPoolManager('conservative');

% Which mRGC/cone mosaic lattice to use
sourceLatticeSizeDegs = 64; 
whichEye = 'right eye';

% Examined spatialChromaticUniformityTradeoff [0: minimize chromatic variance 1: minimize spatial variance]
spatialChromaticUniformityTradeoff = 1.0;

% Which optics to use
whichZernikeDataBase = 'Polans2015';


% First compute-ready mosaic was generated for Polans subject #3
mosaicEccDegs = [0 0]; mosaicSizeDegs = [0.5 0.5];
whichSubjectID = 3;


% Second compute-ready mosaic
mosaicEccDegs = [-3 0]; mosaicSizeDegs = 1*[1 1];

% Extra size (margin)
extraSizeDegs = 1;
mosaicSizeDegs = mosaicSizeDegs+extraSizeDegs*[1 1];

% Use mosaic with RFcenter overlap
employRFCenterOverlappingMosaic = true;

surroundOptimizationStrategy = 'LowerLeftQH1paramsMediumVisualSTFparamTolerance';



% Generate the surroundRetinalConePoolingModel params struct
surroundRetinalConePoolingModelParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(surroundOptimizationStrategy);
if (strcmp(surroundRetinalConePoolingModelParamsStruct.name, 'PackerDacey2002H1FixedCellIndex')) 
    theFixedH1CellIndex = RGCMosaicConstructor.helper.queryUserFor.fixedH1CellIndex(false, []);
    % Add to the surroundRetinalConePoolingModelParamsStruct  the fixed H1 cell index to be used
    surroundRetinalConePoolingModelParamsStruct.fixedH1CellIndex = theFixedH1CellIndex;
end

% Generate the surroundConnectivity simulation params struct
optimizationPositionsAndSizesGridsDummy = [0 0 0 0];
surroundConnectivitySimulationParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundConnectivitySimulationParamsStruct(...
    whichZernikeDataBase, whichSubjectID, employRFCenterOverlappingMosaic, ...
    optimizationPositionsAndSizesGridsDummy, surroundRetinalConePoolingModelParamsStruct);

% Generate filename for the surround-connected, compute-ready MRGCmosaic
surroundConnectedParamsStruct.whichEye = whichEye;
surroundConnectedParamsStruct.eccentricityDegs = mosaicEccDegs;
surroundConnectedParamsStruct.sizeDegs = mosaicSizeDegs;
surroundConnectedParamsStruct.spatialChromaticUniformityTradeoff = spatialChromaticUniformityTradeoff;
surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct = surroundConnectivitySimulationParamsStruct;
surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition = [];
if (employRFCenterOverlappingMosaic)
    surroundConnectivityStage = 'surround connected with center overlap';
else
    surroundConnectivityStage = 'surround connected';
end

% Load the surround-connected, compute-reday MRGCmosaic
theSurroundConnectedMRGCMosaicFullFileName = RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
                surroundConnectedParamsStruct, surroundConnectivityStage);

load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
% Print info on how the surrounds where optimized
%RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.rfSurroundConnectivityParams, 'rfSurroundConnectivityParams')



if (1==2)
% Generate the optics that where used to optimize the mosaic
theOI = theMRGCMosaic.nativeOI(...
    'visualizePSF', true);
end

if (1==1)
% Or, generate alternate optics (adaptive optics)
opticsModification = 'adaptiveOptics6MM';  % To see the valid options type mRGCMosaic.validOpticsModifications
theOI = theMRGCMosaic.nativeOI(...
    'opticsModification', opticsModification, ...
    'visualizePSF', true);
end

if (1==2)
% Or, generate alternate optics, here a custom refraction, -4.0D, appropriate for
% an (x,y) eccentricity of (1, -0.5) degrees
customRefractionDiopters = -4;
opticsModification = 'customRefraction';  % To see the valid options type mRGCMosaic.validOpticsModifications
theOI = theMRGCMosaic.nativeOI(...
    'eccentricityDegs', [0 0], ...
    'opticsModification', opticsModification, ...
    'customRefractionDiopters', customRefractionDiopters, ...
    'visualizePSF', true, ...
    'visualizedWavelengths', 450:20:650);
end


% Compute a demo stimulus scene (checkerboard)
imageFOVdegs = max(theMRGCMosaic.sizeDegs(:));
pixelsPerCheck = 64;
numberOfChecks = 3;
degsPerCheck = imageFOVdegs/numberOfChecks;

theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

% Compute the retinal image
theStimulusRetinalImage = oiCompute(theStimulusScene, theOI);

% Compute the cone mosaic response
[theNoiseFreeConeMosaicExcitationsResponse, theNoisyConeMosaicExcitationsResponses, ~, ~, ...
 theConeMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.inputConeMosaic.compute(...
    theStimulusRetinalImage, ...
    'opticalImagePositionDegs', theMRGCMosaic.eccentricityDegs);

hFig = figure(1);
ax = subplot(1,2,1);
theMRGCMosaic.inputConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'activation', theNoiseFreeConeMosaicExcitationsResponse);

ax = subplot(1,2,2);
theMRGCMosaic.inputConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'activation', theNoisyConeMosaicExcitationsResponses);

[theNoiseFreeMRGCMosaicResponse, ...
 theNoisyMRGCMosaicResponseInstances, ...
 responseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
    theNoiseFreeConeMosaicExcitationsResponse, theConeMosaicResponseTemporalSupportSeconds);



% Set the noise level
theMRGCMosaic.vMembraneGaussianNoiseSigma = std(theNoiseFreeMRGCMosaicResponse(:));

[theNoiseFreeMRGCMosaicResponse, ...
 theNoisyMRGCMosaicResponseInstances, ...
 responseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
    theNoiseFreeConeMosaicExcitationsResponse, theConeMosaicResponseTemporalSupportSeconds);





targetYdegs = theMRGCMosaic.eccentricityDegs(2)+0.25*degsPerCheck;
theROI = regionOfInterest('geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [theMRGCMosaic.eccentricityDegs(1) targetYdegs], ...
            'width', theMRGCMosaic.sizeDegs(1), ...
            'height', 0.03, ...
            'rotation', 0.0...
        ));

hFig = figure(2);
ax = subplot(1,2,1);
theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'activation', theNoiseFreeMRGCMosaicResponse);
hold(ax, 'on');
theROI.visualize('figureHandle', hFig, 'axesHandle', ax);

ax = subplot(1,2,2);
theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'activation', theNoisyMRGCMosaicResponseInstances);
hold(ax, 'on');
theROI.visualize('figureHandle', hFig, 'axesHandle', ax);


visualizedMRGCindices = theROI.indicesOfPointsInside(theMRGCMosaic.rgcRFpositionsDegs);
theVisualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));
[~,idx] = sort(theVisualizedMRGCXcoords, 'ascend');
theVisualizedMRGCXcoords = theVisualizedMRGCXcoords(idx);
visualizedMRGCindices = visualizedMRGCindices(idx);

visualizedConeIndices = theROI.indicesOfPointsInside(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs);
theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
[~,idx] = sort(theVisualizedConeXcoords , 'ascend');
theVisualizedConeXcoords = theVisualizedConeXcoords(idx);
visualizedConeIndices = visualizedConeIndices(idx);

figure(3); clf;
plot(theVisualizedMRGCXcoords, squeeze(theNoiseFreeMRGCMosaicResponse(1,1,visualizedMRGCindices)), 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
hold on
plot(theVisualizedConeXcoords, squeeze(theNoiseFreeConeMosaicExcitationsResponse(1,1,visualizedConeIndices)), 'ko-', 'MarkerSize', 8, 'MarkerFaceColor', [0.75 0.75 0.75]);

