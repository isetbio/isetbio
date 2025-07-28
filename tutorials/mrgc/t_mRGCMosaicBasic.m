
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

    AppleSiliconParPoolManager(4);

    % List the locally-available prebaded mosaics
    % themRGCmosaicFileNames = mRGCMosaic.listPrebakedMosaics();
   
    % Mosaic params
    mosaicParams.eccDegs  = [-4 0];
    mosaicParams.sizeDegs = [3 3];
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = 1;
    mosaicParams.surroundOptimizationSubString = 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';
    mosaicParams.croppedFoVdegs = 1.0;

    % Optics params
    opticsParams.ZernikeDataBase = 'Polans2015';
    opticsParams.subjectRankOrder = 2; 

    % Adaptive optics
    opticsParams.type = 'adaptiveOptics6MM';
    opticsParams.residualWithRespectToNativeOpticsDefocusDiopters = [];

    % native optics + StrehlRatio optimization (what was used to optimize the mosaic
    opticsParams.type = 'nativeOptics';
    opticsParams.residualWithRespectToNativeOpticsDefocusDiopters = [];

    % native optics without StrehlRatio optimization
    %opticsParams.type = 'customRefraction';
    %opticsParams.residualWithRespectToNativeOpticsDefocusDiopters = 0.0;
    opticsParams.visualizePSFonTopOfConeMosaic = true;

    % Load the desired mRGCmosaic and the associated optics
    [theMRGCmosaic, theOI] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);
    theMRGCmosaic.visualize();

    % Input stimulus
    imageFOVdegs = min(theMRGCmosaic.sizeDegs);
    pixelsPerCheck = 256;
    numberOfChecks = 5;

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


