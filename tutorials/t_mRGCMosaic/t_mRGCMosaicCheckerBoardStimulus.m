function t_mRGCMosaicCheckerBoardStimulus

    % Specify the precomputed mosaic's eccentricity
    horizontalEccDegs = 7;

    % Load the precomputed mRGCMosaic
    theMRGCMosaic = loadPreComputedMRGCMosaic(horizontalEccDegs);

    % Noise-free responses
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'none';
    theMRGCMosaic.noiseFlag = 'none';

    % Input stimulus
    imageFOVdegs = 7.0;
    pixelsPerCheck = 64;
    numberOfChecks = 4;

    % Compute the stimulus scene
    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

    % Compute the retinal image
    theOI = theMRGCMosaic.theNativeOptics;
    theStimulusRetinalImage = oiCompute(theStimulusScene, theOI);

    % Compute the cone mosaic response
    theConeMosaicResponse = theMRGCMosaic.inputConeMosaic.compute(...
        theStimulusRetinalImage, ...
        'opticalImagePositionDegs', [7. 0]);

    % Compute the mRGC mosaic response using the default RGCRGgains
    theConeMosaicResponseTemporalSupportSeconds = [0];
    theMRGCMosaicResponseDefaultGain = theMRGCMosaic.compute( ...
             theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);


    % Retrieve the default rgcRFgains
    defaultRGCRFgains = theMRGCMosaic.rgcRFgains;
    max(theMRGCMosaic.rgcRFgains)
    min(theMRGCMosaic.rgcRFgains)
    pause

    % Set the rgcRFgains as 1 / center-integrated retinal cone apertures
    method = '1/integrated center retinal cone apertures';
    maxGain = 8.5e-12;
    methodParams = maxGain;
    theMRGCMosaic.setPeakGains(method, methodParams);

    % Compute the mRGC mosaic response using the integrated cone aperture gain
    theMRGCMosaicResponseIntegratedConeApertureBasedGain = theMRGCMosaic.compute( ...
             theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);

    % Compute indices of visualized cones and RGCs along the Y = -0.65 degs
    targetYdegs = -0.7;
    
    [visualizedConeIndices, theVisualizedConeXcoords, ...
     visualizedMRGCindices, theVisualizedMRGCXcoords, theROI] = extractVisualizedConeAndRGCindices(theMRGCMosaic, targetYdegs);

    % Extract responses of visualized cones
    theVisualizedConeResponses = squeeze(theConeMosaicResponse(1,1,visualizedConeIndices));

    % Extract responses of visualized mRGCs
    theVisualizedMRGCResponsesDefaultGain = squeeze(theMRGCMosaicResponseDefaultGain(1,1,visualizedMRGCindices));
    theVisualizedMRGCResponsesIntegratedConeApertureBasedGain = squeeze(theMRGCMosaicResponseIntegratedConeApertureBasedGain(1,1,visualizedMRGCindices));

    % Restore the default rgcRFgains
    method = 'arbitrary';
    methodParams = defaultRGCRFgains;
    theMRGCMosaic.setPeakGains(method, methodParams);


    % Visualize the mRGCMosaic
    visualizeTheMRGCMosaic(1,theMRGCMosaic);

    % Visualize the activation of the input cone mosaic
    visualizeRetinalOpticalImage(2, theMRGCMosaic, theStimulusRetinalImage);

    % Visualize the activation of the input cone mosaic
    visualizedResponseGain = 0.8;

    visualizeInputConeMosaicActivation(3, theMRGCMosaic, theConeMosaicResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        visualizedResponseGain, numberOfChecks);

    % Visualize the activation of the mRGCmosaic (default gain)
    visualizeMRGCMosaicActivation(4, theMRGCMosaic, theMRGCMosaicResponseDefaultGain, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponsesDefaultGain, ...
        visualizedResponseGain, numberOfChecks, 'defaultGain');

    % Visualize the activation of the mRGCmosaic (integrated cone apertures based gain)
    visualizeMRGCMosaicActivation(5, theMRGCMosaic, theMRGCMosaicResponseIntegratedConeApertureBasedGain , ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponsesIntegratedConeApertureBasedGain, ...
        visualizedResponseGain, numberOfChecks, 'integratedConeAperturesBasedGain');
end




% HELPER FUNCTIONS

function visualizeMRGCMosaicActivation(figNo, theMRGCMosaic, theMRGCMosaicResponse, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponses, ...
        visualizedResponseGain, numberOfChecks, postFix)

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [0 0 0], 'Name', postFix);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);

    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theMRGCMosaicResponse, ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0 0 0], ...
        'domainVisualizationLimits', [theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1] theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1]], ...
        'domainVisualizationTicks', struct('x', 4:10, 'y', -5:5));

    hold(ax, 'on');
    if (1==1)
        theROI.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'xLims', theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1], ...
            'yLims', theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1], ...
            'fillColor', [0.8 0.8 0.8 0.4]);

        plot(ax, theVisualizedMRGCXcoords, targetYdegs + theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:))*visualizedResponseGain, 'r-', 'LineWidth', 1.5);
        plot(ax, theVisualizedMRGCXcoords, targetYdegs + theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:))*visualizedResponseGain, 'ro', 'MarkerFaceColor', [1 0.8 0.5], 'LineWidth', 1.5);
        plot(ax, theVisualizedMRGCXcoords, targetYdegs + theVisualizedMRGCResponses/max(theMRGCMosaicResponse(:))*0.0, 'r-', 'LineWidth', 1.5);
    end

    set(hFig, 'Color', [0 0 0]);
    set(ax, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', 30);
    grid(ax, 'off')
    NicePlot.exportFigToPDF(sprintf('theMRGCMosaicActivation_%dChecks_%s.pdf', numberOfChecks, postFix), hFig, 300);
end


function visualizeInputConeMosaicActivation(figNo, theMRGCMosaic, theConeMosaicResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        visualizedResponseGain, numberOfChecks)

    hFig = figure(figNo);clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [0 0 0]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);

    theMRGCMosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theConeMosaicResponse, ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0 0 0], ...
        'domainVisualizationLimits', [theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1] theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1]], ...
        'domainVisualizationTicks', struct('x', 4:10, 'y', -5:5));

    hold(ax, 'on');
    theROI.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'xLims', theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1], ...
        'yLims', theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1], ...
        'fillColor', [0.8 0.8 0.8 0.4]);

    plot(ax, theVisualizedConeXcoords, targetYdegs + theVisualizedConeResponses/max(theConeMosaicResponse(:))*visualizedResponseGain, 'r-', 'LineWidth', 1.5);
    plot(ax, theVisualizedConeXcoords, targetYdegs + theVisualizedConeResponses/max(theConeMosaicResponse(:))*visualizedResponseGain, 'ro', 'MarkerFaceColor', [1 0.8 0.5], 'LineWidth', 1.5);
    plot(ax, theVisualizedConeXcoords, targetYdegs + theVisualizedConeResponses/max(theConeMosaicResponse(:))*0.0, 'r-', 'LineWidth', 1.5);
    
    set(hFig, 'Color', [0 0 0]);
    set(ax, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', 30);
    grid(ax, 'off')
    NicePlot.exportFigToPDF(sprintf('theConeMosaicActivation_%dChecks.pdf', numberOfChecks), hFig, 300);
end


function visualizeRetinalOpticalImage(figNo, theMRGCMosaic, theStimulusRetinalImage)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);

    theMRGCMosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'conesAlpha', 0.0, ...
        'conesEdgeAlpha', 0.5, ...
        'withSuperimposedOpticalImage', theStimulusRetinalImage, ...
        'superimposedOIAlpha', 0.7, ...
        'domainVisualizationLimits', [theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1] theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1]], ...
        'domainVisualizationTicks', struct('x', 4:10, 'y', -5:5));
    NicePlot.exportFigToPDF('theRetinalOpticalImage.pdf', hFig, 300);

end


function visualizeTheMRGCMosaic(figNo,theMRGCMosaic)

    % Here, we visualize the mRGCMosaic (RF centers) without with the input cone mosaic.
    % In this visualization, the gray contours identify the mRGC RF centers
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyInputCones', ~true, ...
        'identifyPooledCones', ~true, ...
        'backgroundColor', [1 1 1], ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'domainVisualizationLimits', [theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1] theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1]], ...
         'domainVisualizationTicks', struct('x', 4:10, 'y', -5:5), ...
         'plotTitle', ' ');
    set(hFig, 'Color', [1 1 1]);
    set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'LineWidth', 1.5, 'FontSize', 30);
    grid(ax, 'off')
    NicePlot.exportFigToPDF('theMRGCMosaic.pdf', hFig, 300)
end



function [visualizedConeIndices, theVisualizedConeXcoords, ...
          visualizedMRGCindices, theVisualizedMRGCXcoords, theROI] = extractVisualizedConeAndRGCindices(theMRGCMosaic, targetYdegs)

    % Define an ROI
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [theMRGCMosaic.eccentricityDegs(1) targetYdegs], ...
            'width', theMRGCMosaic.sizeDegs(1), ...
            'height', 0.2, ...
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

function theMRGCMosaic = loadPreComputedMRGCMosaic(horizontalEccDegs)

    % Choose one of the computed mRGCMosaics
    % 1. mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(horizontalEccDegs);

    % 2. optics params
    opticsParams = MosaicPoolingOptimizer.getOpticsParams(mosaicParams);

    % 3. retinal RF pooling params
    retinalRFmodelParams = MosaicPoolingOptimizer.getSurroundParams(mosaicParams, opticsParams);

    theMRGCMosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
        mosaicParams, opticsParams, retinalRFmodelParams);
end
