
%% Introduction to the 2025 midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates the following:
%    (i)   how to load apre-baked midget RGC mosaic, 
%    (ii)  generate custom optics for the computation,
%    (iii) compute the response of the mRGC mosaic to a checkberboard stimulus, and
%    (iv)  visualize the responses of the mRGC mosaic and of its input cone mosaic
%


% History:
%    07/28/25  NPC  Wrote it.

function t_mRGCMosaicBasic
    %% Close all figures
    close all;

    AppleSiliconParPoolManager(4);

    % List the locally-available prebaded mosaics
    % themRGCmosaicFileNames = mRGCMosaic.listPrebakedMosaics();
   
    % Load one of the prebaked mRGC mosaics. 
    % We need to specify 3 arguments:
    % (i) the eccentricity and size of the mosaic
    % (ii) the surround optimization method
    % (iii) the optics under which the mosaic was optimized
    % (iv) the optics to employ for the computation at hand

    % (i) Eccentricity (-4.0) and size (3x3deg) 
    mosaicParams.eccDegs  = [-4 0];
    mosaicParams.sizeDegs = [3 3];

    % (ii) Surround optimization method
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = 1;
    mosaicParams.surroundOptimizationSubString = 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';

    % We will employ a smaller (2x1) patch centered at (-5,0), cropped from
    % the full mosaic. If cropParams is not specified, we employ the full
    % mosaic
    mosaicParams.cropParams = struct(...
        'sizeDegs', [2 1], ...
        'eccentricityDegs', [-5 0] ...
        );


    % (iii) Optics under which the mosaic was optimized
    opticsParams.ZernikeDataBase = 'Polans2015';
    opticsParams.subjectRankOrder = 2; 
    opticsParams.visualizePSFonTopOfConeMosaic = true;


    % (iv) Optics to employ for the computation at hand
    
    % EITHER adaptive optics (diffraction limited with 6 mm pupil)
    %opticsParams.type = 'adaptiveOptics6MM';
    %opticsParams. refractiveErrorDiopters = [];
    
    % OR the native optics + StrehlRatio optimization (what was used to optimize the mosaic
    %opticsParams.type = 'nativeOptics';
    %opticsParams.refractiveErrorDiopters = [];

    % OR the native optics without a custom refraction, here -3.5D
    opticsParams.type = 'customRefraction';
    opticsParams.refractiveErrorDiopters = -3.5;
    
    
    % Load the desired mRGCmosaic and generated the optics for the
    % computation
    [theMRGCmosaic, theOpticsToEmploy] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);
    theMRGCmosaic.visualize();

    % Input stimulus
    imageFOVdegs = min(theMRGCmosaic.sizeDegs);
    pixelsPerCheck = 256;
    numberOfChecks = 5;

    % Compute the stimulus scene
    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

    % Compute the retinal image
    theStimulusRetinalImage = oiCompute(theOpticsToEmploy, theStimulusScene);

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
