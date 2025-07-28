
%% Introduction to the 2025 midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates: how to load of of pre-baked midget RGC mosaics, compute
%    its response to a checkberboard stimulus, and visualize the response
%    of the input cone mosaic along with the response of the mRGC mosaic
%


% History:
%    07/28/25  NPC  Wrote it.

function t_mRGCMosaicBasic
    %% Close all figures
    close all;

    % Load an mRGCmosaic located the far periphery
    visualizePSFonTopOfConeMosaic = true;
    [theMRGCmosaic, theOI] = farPeripheryMRGCmosaiAndOptics(visualizePSFonTopOfConeMosaic);
    theMRGCmosaic.visualize();

    % Input stimulus
    imageFOVdegs = min(theMRGCmosaic.sizeDegs);
    pixelsPerCheck = 256;
    numberOfChecks = 4;

    % Compute the stimulus scene
    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

    % Compute the retinal image
    theStimulusRetinalImage = oiCompute(theOI, theStimulusScene);

    % Compute the input cone mosaic response to the retinal image of the stimulus
    [theNoiseFreeConeMosaicExcitationsResponse, ...
     theNoisyConeMosaicExcitationsResponses, ~, ~, ...
     theConeMosaicResponseTemporalSupportSeconds] = theMRGCmosaic.inputConeMosaic.compute(...
        theStimulusRetinalImage, ...
        'opticalImagePositionDegs', theMRGCmosaic.eccentricityDegs);

    % Compute the mRGCmosaic response to the input cone mosaic response
    mRGCNonLinearityParams = [];
    [theNoiseFreeSpatioTemporalMRCMosaicResponse, ~, ...
     theMRGCMosaicResponseTemporalSupportSeconds] = theMRGCmosaic.compute( ...
                theNoiseFreeConeMosaicExcitationsResponse, ...
                theConeMosaicResponseTemporalSupportSeconds, ...
                'nonLinearityParams', mRGCNonLinearityParams);

    hFig = figure(1);
    set(hFig, 'Position', [10 10 1980 910]);

    ax = subplot(1,2,1);
    theMRGCmosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theNoiseFreeConeMosaicExcitationsResponse, ...
        'plotTitle', 'input cone mosaic response');

    ax = subplot(1,2,2);
    theMRGCmosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theNoiseFreeSpatioTemporalMRCMosaicResponse, ...
        'plotTitle', 'mRGC mosaic response');


end

% Supporting functions
function [theMRGCMosaic, theOI] = farPeripheryMRGCmosaiAndOptics(visualizePSFonTopOfConeMosaic)

    theOpticsSubject = 'Polans2015-2';
    theMosaicXYeccentricityDegs = [-32.0 0.0];
    theMosaicXYsizeDegs = [9 9];

    prebakedMRGCMosaicDir = 'isettools/ganglioncells/data/prebakedRGCmosaics/ONmRGCmosaics';
    spatialCompactnessSpectralPurityTradeoff = 1;
    opticsSubString = sprintf('Optics_%s_maxStrehlRatio', theOpticsSubject);
    surroundOptimizationSubString = 'PackerDacey2002H1freeUpperH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';

    mRGCMosaicFilename = sprintf('MRGCMosaic_RE_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_Phi_%1.2f_%s_srndModel_%s.mat', ...
        theMosaicXYeccentricityDegs(1), theMosaicXYeccentricityDegs(2), ...
        theMosaicXYsizeDegs(1), theMosaicXYsizeDegs(2), ...
        spatialCompactnessSpectralPurityTradeoff, opticsSubString, surroundOptimizationSubString);

    load(fullfile(isetbioRootPath, prebakedMRGCMosaicDir,mRGCMosaicFilename), 'theMRGCMosaic');

    % Employ the native optics (what was used to optimize the surround)
    opticsForSTFresponses = 'nativeOptics';  % native optics + StrehlRatio optimization (what was used to optimize the mosaic)
    
    %opticsForSTFresponses = 'adaptiveOptics6MM';

    opticsForSTFresponses = 'customRefraction';  % native optics without StrehlRatio optimization )
    residualWithRespectToNativeOpticsDefocusDiopters = 0.0;

    

    % Generate the optics for the mosaic
    [theOI, thePSF] = RGCMosaicAnalyzer.compute.opticsForResponses(...
        theMRGCMosaic, opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, visualizePSFonTopOfConeMosaic);
 end


