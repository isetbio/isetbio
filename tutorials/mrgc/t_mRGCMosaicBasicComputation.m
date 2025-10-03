%% Introduction to the 2025 midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates the following:
%    (i)   how to load a pre-baked midget RGC mosaic, 
%    (ii)  generate different optics for the computation,
%    (iii) compute the response of the mRGC mosaic to a checkberboard stimulus, and
%    (iv)  visualize the responses of the mRGC mosaic and of its input cone mosaic
%


% History:
%    07/28/25  NPC  Wrote it.

function t_mRGCMosaicBasicComputation
    % Initialize
    clearvars;
    close all;

    % Use 6 CPU cores
    AppleSiliconParPoolManager(4);

    % List the locally-available prebaded mosaics
    % themRGCmosaicFileNames = mRGCMosaic.listPrebakedMosaics();
   
    % Load one of the prebaked mRGC mosaics. 
    % We need to specify 4 pieces of information:
    % (A) the eccentricity and size of the mosaic
    % (B) the surround optimization method
    % (C) the optics under which the mosaic was optimized

    % (A) Eccentricity
    % Load an 11x11 deg mosaic that was synthesized at 7 degrees along the nasal meridian
    % using optics from Polans subject with rank order 9
    mosaicParams.eccDegs  = [7 0];
    mosaicParams.sizeDegs = [11 11];

    
    % (B) Surround optimization method
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = 1;
    mosaicParams.surroundOptimizationSubString = 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';


    % (C) Optics under which the mosaic was optimized
    opticsParams.ZernikeDataBase = 'Polans2015';
    opticsParams.subjectRankOrder = 9; 
    opticsParams.type = 'nativeOptics';
    opticsParams.refractiveErrorDiopters = [];
    opticsParams.visualizePSFonTopOfConeMosaic = ~true;


    % (D) Crop a patch (8x4 degs) of the mosaic, centered at 7 degrees
    mosaicParams.cropParams = struct(...
        'sizeDegs', [6 3], ...
        'eccentricityDegs', [7 0] ...
        );

    % Append to opticsParams, information on which optics to employ for the computation at hand
    
    % EITHER adaptive optics (diffraction limited with 6 mm pupil)
    %opticsParams.type = 'adaptiveOptics6MM';
    %opticsParams. refractiveErrorDiopters = [];
    
    % OR the native optics + StrehlRatio optimization (what was used to optimize the mosaic
    opticsParams.type = 'nativeOptics';
    opticsParams.refractiveErrorDiopters = [];

    % OR the native optics without a custom refraction, here -3.5D
    %opticsParams.type = 'customRefraction';
    %opticsParams.refractiveErrorDiopters = -3.5;
    
    % OR the native optics without a custom refraction, here -0.5D
    %opticsParams.type = 'customRefraction';
    %opticsParams.refractiveErrorDiopters = -0.5;
    
    % Load the desired mRGCmosaic and generated the optics for the computation
    [theMRGCmosaic, theOpticsToEmploy, thePSF] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);

    % Subdirectory for exporting the generated PDFs
    exportVisualizationPDFdirectory = 'mosaicComputePDFs';
    
    % The full path where the generated PDFs will be stored can be found like so:
    % p = getpref('isetbio');
    % fullPath = fullfile(p.rgcResources.figurePDFsDir,exportVisualizationPDFdirectory);

    % Visualize the mosaic of mRGC RF centers
    % identifying cones that are pooled by the RF center mechanism with
    % a weight >= mRGCMosaic.sensitivityAtPointOfOverlap;
    % This representation is like the representation used in visualizing 
    % mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
    minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;
    
    % Visualization limits and ticks
    visualizedWidthDegs = theMRGCmosaic.inputConeMosaic.sizeDegs(1);
    visualizedHeightDegs = theMRGCmosaic.inputConeMosaic.sizeDegs(2);;
    domainVisualizationLimits(1:2) = theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 1];
    domainVisualizationLimits(3:4) = theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1];
    domainVisualizationTicks = struct(...
        'x', theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 -0.5 0 0.5 1], ...
        'y', theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedWidthDegs * [-1 -0.5 0 0.5 1]);

    ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');

    hFig = figure(1); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    theMRGCmosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'inputConesAlpha', 0.5, ...
        'identifiedConeAperture', 'lightCollectingArea4sigma', ...
        'identifiedConeApertureThetaSamples', 16, ...
        'minConeWeightVisualized', minCenterConeWeight, ...
        'centerSubregionContourSamples', 32, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
        'visualizationPDFfileName', sprintf('mRGCmosaicMinCenterConeWeight_%2.3f', minCenterConeWeight), ...
        'exportVisualizationPDF', true, ...
        'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);


    % Input stimulus: match its size to the input cone mosaic
    imageFOVdegs = max(theMRGCmosaic.inputConeMosaic.sizeDegs);
    pixelsPerCheck = 128;
    numberOfChecks = 5;

    % Compute the stimulus scene
    fprintf('Generating checkerboard scene. Please wait ...\n')
    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);


    % Compute the retinal image
    fprintf('Computing retinal image. Please wait ...\n')
    theStimulusRetinalImage = oiCompute(theOpticsToEmploy, theStimulusScene);

    % Empty some RAM
    clear 'theStimulusScene';

    % Compute the input cone mosaic response to the retinal image of the stimulus
    fprintf('Computing input cone mosaic response. Please wait ...\n')
    [theNoiseFreeConeMosaicExcitationsResponse, ...
     theNoisyConeMosaicExcitationsResponses, ~, ~, ...
     theConeMosaicResponseTemporalSupportSeconds] = theMRGCmosaic.inputConeMosaic.compute(...
        theStimulusRetinalImage, ...
        'opticalImagePositionDegs', theMRGCmosaic.eccentricityDegs);

    % Empty some RAM
    clear 'theStimulusRetinalImage';


    % Compute the mRGCmosaic response to the input cone mosaic response
    fprintf('Computing mRGC mosaic response. Please wait ...\n')
    [theNoiseFreeSpatioTemporalMRCMosaicResponse, ~, ...
     theMRGCMosaicResponseTemporalSupportSeconds] = theMRGCmosaic.compute( ...
                theNoiseFreeConeMosaicExcitationsResponse, ...
                theConeMosaicResponseTemporalSupportSeconds);


    % Visualize the response of the input cone mosaic
    ff = PublicationReadyPlotLib.figureComponents('1x2 standard figure');
    ff.backgroundColor = [0 0 0];
    hFig = figure(2); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax1 = theAxes{1,1};
    ax2 = theAxes{1,2};

    theMRGCmosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax1, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'activation', theNoiseFreeConeMosaicExcitationsResponse, ...
        'plotTitle', 'input cone mosaic response');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax1,ff);


    % Visualize the response of the mRGCmosaic
    theMRGCmosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax2, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'activation', theNoiseFreeSpatioTemporalMRCMosaicResponse, ...
        'plotTitle', 'mRGC mosaic response');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax2,ff);


    % Compute indices of visualized cones and RGCs along the Y = 0.15 degs
    targetYdegs = -0.42;
    targetHeightDegs = 0.05;

    % Compute indices of mRGCs and input cones that are in the target ROI
    [visualizedConeIndices, theVisualizedConeXcoords, ...
     visualizedMRGCindices, theVisualizedMRGCXcoords, theROI] = extractVisualizedConeAndRGCindices(theMRGCmosaic, targetYdegs, targetHeightDegs);


    % Extract responses of visualized cones
    theVisualizedConeResponses = squeeze(theNoiseFreeConeMosaicExcitationsResponse(1,1,visualizedConeIndices));

    % Extract responses of visualized mRGCs
    theVisualizedMRGCResponses = squeeze(theNoiseFreeSpatioTemporalMRCMosaicResponse(1,1,visualizedMRGCindices));
    

    % Visualize the mRGC mosaic activation
    ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');
    ff.backgroundColor = [0 0 0];

    hFig = figure(3); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    visualizeMRGCMosaicActivationWithinROI(hFig, ax, theMRGCmosaic, theNoiseFreeSpatioTemporalMRCMosaicResponse, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponses, ...
        domainVisualizationLimits, domainVisualizationTicks);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export PDF
    p = getpref('isetbio');
    figureDir = fullfile(p.rgcResources.figurePDFsDir, exportVisualizationPDFdirectory);
    NicePlot.exportFigToPDF(fullfile(figureDir,'checkerboardMRGCMosaicActivation.pdf'), hFig, 300);


    % Visualize the input cone mosaic activation
    hFig = figure(4); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    visualizeInputConeMosaicActivation(hFig, ax, theMRGCmosaic, theNoiseFreeConeMosaicExcitationsResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        domainVisualizationLimits, domainVisualizationTicks);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export PDF
    p = getpref('isetbio');
    figureDir = fullfile(p.rgcResources.figurePDFsDir, exportVisualizationPDFdirectory);
    NicePlot.exportFigToPDF(fullfile(figureDir,'checkerboardInputConeMosaicActivation.pdf'), hFig, 300);


end


%
% HELPER FUNCTIONS
%



function visualizeMRGCMosaicActivationWithinROI(hFig, ax, theMRGCMosaic, theMRGCMosaicResponse, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponses, ...
        domainVisualizationLimits, domainVisualizationTicks)


    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theMRGCMosaicResponse, ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0 0 0], ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks);

    hold(ax, 'on');
    
    theROI.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'xLims', theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1], ...
        'yLims', theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1], ...
        'fillColor', [0.8 0.8 0.8 0.4]);


    plot(ax, theVisualizedMRGCXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:)), 'k-', 'LineWidth', 5);
    plot(ax, theVisualizedMRGCXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:)), 'r-', 'LineWidth', 4);
    plot(ax, theVisualizedMRGCXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:)), 'ro', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.8 0.5], 'LineWidth', 1.5);
    plot(ax, theVisualizedMRGCXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:))*0.0, 'r-', 'LineWidth', 2);


    set(ax, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', 30);
    set(ax, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));
    set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
    grid(ax, 'off')
end


function visualizeInputConeMosaicActivation(hFig, ax, theMRGCMosaic, theConeMosaicResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        domainVisualizationLimits, domainVisualizationTicks)


    theMRGCMosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theConeMosaicResponse, ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0 0 0], ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks);

    hold(ax, 'on');
    theROI.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'xLims', theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1], ...
        'yLims', theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1], ...
        'fillColor', [0.8 0.8 0.8 0.4]);

    plot(ax, theVisualizedConeXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedConeResponses/max(theConeMosaicResponse(:)), 'k-', 'LineWidth', 5);
    plot(ax, theVisualizedConeXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedConeResponses/max(theConeMosaicResponse(:)), 'r-', 'LineWidth', 4);
    plot(ax, theVisualizedConeXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedConeResponses/max(theConeMosaicResponse(:)), 'ro', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.8 0.5], 'LineWidth', 1.5);
    plot(ax, theVisualizedConeXcoords, targetYdegs + theMRGCMosaic.sizeDegs(2)*0.4 * theVisualizedConeResponses/max(theConeMosaicResponse(:))*0.0, 'r-', 'LineWidth', 2);
    
    set(hFig, 'Color', [0 0 0]);
    set(ax, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', 30);
    set(ax, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));
    set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
    grid(ax, 'off')

end



function [visualizedConeIndices, theVisualizedConeXcoords, ...
          visualizedMRGCindices, theVisualizedMRGCXcoords, theROI] = extractVisualizedConeAndRGCindices(theMRGCMosaic, targetYdegs, targetHeightDegs)

    % Define an ROI
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [theMRGCMosaic.eccentricityDegs(1) targetYdegs], ...
            'width', theMRGCMosaic.sizeDegs(1), ...
            'height', targetHeightDegs, ...
            'rotation', 0.0...
        ));

    visualizedConeIndices = theROI.indicesOfPointsInside(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs);
    visualizedMRGCindices = theROI.indicesOfPointsInside(theMRGCMosaic.rgcRFpositionsDegs);

    theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
    theVisualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));

    [~,idx] = sort(theVisualizedConeXcoords, 'ascend');
    visualizedConeIndices = visualizedConeIndices(idx);
    % Exclude S-cones
    idx = find(theMRGCMosaic.inputConeMosaic.coneTypes(visualizedConeIndices) == cMosaic.SCONE_ID);
    [~, idx] = setdiff(visualizedConeIndices, visualizedConeIndices(idx));
    visualizedConeIndices = visualizedConeIndices(idx);
    theVisualizedConeXcoords = squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));

    [~,idx] = sort(theVisualizedMRGCXcoords, 'ascend');
    theVisualizedMRGCXcoords = theVisualizedMRGCXcoords(idx);
    visualizedMRGCindices = visualizedMRGCindices(idx);
end

