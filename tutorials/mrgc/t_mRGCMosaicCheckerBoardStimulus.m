function t_mRGCMosaicCheckerBoardStimulus(options)
% Compute the response of an mRGCmosaic to a checkerboard stimulus
%
% Requires locally configured precomputed mRGC resources.
% SkipFile
%
% History:    
%    7/26/2026 NPC  Modernized it

% Examples:
%{

    t_mRGCMosaicCheckerBoardStimulus();

    t_mRGCMosaicCheckerBoardStimulus(...
        'rgcMosaicName', 'JCNpaperTemporal7DegsMosaic', ...
        'opticsSubjectName', 'JCNpaperDefaultSubject');

%}

arguments

    % ---- Mosaic specifiers for selecting a prebaked mRGC mosaic ---

    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'JCNpaperNasal2DegsTinyMosaic';


    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.coneMosaicSpecies  (1,:) char {mustBeMember(options.coneMosaicSpecies,{'human','macaque'})} = 'human';

     % ----- Which subject optics to employ -----
    options.opticsSubjectName (1,:) ...
        char ...
        {...
        mustBeMember(options.opticsSubjectName, ...
            { ...
            'JCNpaperDefaultSubject' ...
            'JCNpaperSecondSubject' ...
            'VSS2024TalkFirstSubject' ...
            'VSS2024TalkSecondSubject' ...
            'JCNpaperStrehlRatio_0.87' ...
            'JCNpaperStrehlRatio_0.72' ...
            'JCNpaperStrehlRatio_0.59' ...
            'JCNpaperStrehlRatio_0.60' ...
            'JCNpaperStrehlRatio_0.27' ...
            'JCNpaperStrehlRatio_0.23' ...
            'JCNpaperStrehlRatio_0.21' ...
            'JCNpaperStrehlRatio_0.19' ...
            'JCNpaperStrehlRatio_0.09' ...
            } ...
            ) ...
        } ...
        = 'JCNpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % See RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct
    % for all existing options
    options.targetVisualSTFdescriptor (1,:) char = 'default';


    % ------ Visualization options ----
    % Visualize cone pooling maps for a target RGC
    options.targetRGCindexForVisualizingConePoolingMaps (1,:) double = [];

    % Whether to generate a video of RFpooling maps along the horizontal meridian
    options.generateVideoOfRFpoolingMapsAlongHorizontalMeridian (1,1) logical = false;

    % Whether to close previously open figures
    options.closePreviouslyOpenFigures (1,1) logical = true;
end


    % Set flags from key/value pairs
    
    % Mosaic specifiers for selecting a prebaked mRGC mosaic
    rgcMosaicName = options.rgcMosaicName;
    coneMosaicSpecies = options.coneMosaicSpecies;
    opticsSubjectName = options.opticsSubjectName;
    targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;



    %% Close all figures
    close all;

    % Configure a conservative parpool manager. This gives at least 8 GB RAM/core
    % You might need to set this to true if you are hitting a strange
    % parpool crash.  Setting it to true, however, invokes a new parpool
    % open and restore which slows things down in most cases.
    fancyParpoolControl = false;
    if (fancyParpoolControl)
        ASPPManager = AppleSiliconParPoolManager('conservative');
    end

    % Control saving of figures.  We don't want tutorials
    % saving things into the isetbio source tree.
    saveFigures = true;
    figureDir = figExporter.figureDir(mfilename, saveFigures);

    % Load the mRGCMosaic and the optics that were used to optimize it
    % See t_mRGCMosaicBasic for why 'wavefrontSpatialSamples' is set below
    % the default of 501.
    [theMRGCMosaic, theOI] = mRGCMosaic.loadPrebakedMosaic(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
            'computeTheMosaicOptics', true, ...
            'wavefrontSpatialSamples', 301);
   

    % Noise-free responses
    theMRGCMosaic.inputConeMosaic.noiseFlag = 'none';
    theMRGCMosaic.noiseFlag = 'none';

    % Input stimulus
    imageFOVdegs = 7.0;
    % 128 pixels per check gives a 1024x1024 scene, i.e. ~146 pixels/deg at
    % the 7 deg FOV set below. Raising this to 256 quadruples the scene and
    % costs about 8 seconds more in oiCompute and cone excitation, while
    % changing the mRGC responses by ~0.3%.
    pixelsPerCheck = 128;
    numberOfChecks = 4;

    % Compute the stimulus scene
    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

    % Compute the retinal image
    theStimulusRetinalImage = oiCompute(theOI,theStimulusScene);

    % Compute the cone mosaic response
    theConeMosaicResponse = theMRGCMosaic.inputConeMosaic.compute(...
        theStimulusRetinalImage, ...
        'opticalImagePositionDegs', theMRGCMosaic.eccentricityDegs);

    % Compute the mRGC mosaic response using the default RGCRGgains
    theConeMosaicResponseTemporalSupportSeconds = [0];
    theMRGCMosaicResponseDefaultGain = theMRGCMosaic.compute( ...
             theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);


    % Retrieve the default rgcRFgains
    defaultRGCRFgains = theMRGCMosaic.responseGains;

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
    visualizeTheMRGCMosaic(1,theMRGCMosaic, saveFigures, figureDir);

    % Visualize the activation of the input cone mosaic
    visualizeRetinalOpticalImage(2, theMRGCMosaic, theStimulusRetinalImage, saveFigures, figureDir);

    % Visualize the activation of the input cone mosaic
    visualizedResponseGain = 0.8;

    visualizeInputConeMosaicActivation(3, theMRGCMosaic, theConeMosaicResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        visualizedResponseGain, numberOfChecks, saveFigures, figureDir);

    % Visualize the activation of the mRGCmosaic (default gain)
    visualizeMRGCMosaicActivation(4, theMRGCMosaic, theMRGCMosaicResponseDefaultGain, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponsesDefaultGain, ...
        visualizedResponseGain, numberOfChecks, 'defaultGain', saveFigures, figureDir);

    % Visualize the activation of the mRGCmosaic (integrated cone apertures based gain)
    visualizeMRGCMosaicActivation(5, theMRGCMosaic, theMRGCMosaicResponseIntegratedConeApertureBasedGain , ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponsesIntegratedConeApertureBasedGain, ...
        visualizedResponseGain, numberOfChecks, 'integratedConeAperturesBasedGain', saveFigures, figureDir);

    % Restore previous parpool config
    if (fancyParpoolControl)
        ASPPManager.restoreLastParpoolSize();
    end
end

% HELPER FUNCTIONS
function visualizeMRGCMosaicActivation(figNo, theMRGCMosaic, theMRGCMosaicResponse, ...
        theROI, targetYdegs, theVisualizedMRGCXcoords, theVisualizedMRGCResponses, ...
        visualizedResponseGain, numberOfChecks, postFix, saveFigures, figureDir)

    hFig = figure(figNo); clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [0 0 0], 'Name', postFix);
    hFig.Units = originalFigureUnits;
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
    if (saveFigures)
         NicePlot.exportFigToPDF(fullfile(figureDir,sprintf('theMRGCMosaicActivation_%dChecks_%s.pdf', numberOfChecks, postFix)), hFig, 300);
    end
end


function visualizeInputConeMosaicActivation(figNo, theMRGCMosaic, theConeMosaicResponse, ...
        theROI, targetYdegs, theVisualizedConeXcoords, theVisualizedConeResponses, ...
        visualizedResponseGain, numberOfChecks, saveFigures, figureDir)

    hFig = figure(figNo);clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [0 0 0]);
    hFig.Units = originalFigureUnits;
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
    if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(figureDir,sprintf('theConeMosaicActivation_%dChecks.pdf', numberOfChecks)), hFig, 300);
    end
end


function visualizeRetinalOpticalImage(figNo, theMRGCMosaic, theStimulusRetinalImage, saveFigures, figureDir)
    hFig = figure(figNo); clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    hFig.Units = originalFigureUnits;
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);

    theMRGCMosaic.inputConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'conesAlpha', 0.0, ...
        'conesEdgeAlpha', 0.5, ...
        'outlinedConesWithIndices', 1:theMRGCMosaic.inputConeMosaic.conesNum, ...
        'withSuperimposedOpticalImage', theStimulusRetinalImage, ...
        'withsuperimposedOpticalImageAlpha', 0.7, ...
        'domainVisualizationLimits', [theMRGCMosaic.eccentricityDegs(1) + theMRGCMosaic.sizeDegs(1)*0.51*[-1 1] theMRGCMosaic.eccentricityDegs(2) + theMRGCMosaic.sizeDegs(2)*0.51*[-1 1]], ...
        'domainVisualizationTicks', struct('x', 4:10, 'y', -5:5));
    if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(figureDir,'theRetinalOpticalImage.pdf'), hFig, 300);
    end

    exportDir = fullfile(isetbioRootPath,'local');
    saveas(hFig,fullfile(exportDir,'test.tiff'),'tiff');
end


function visualizeTheMRGCMosaic(figNo,theMRGCMosaic, saveFigures, figureDir)

    % Here, we visualize the mRGCMosaic (RF centers) without with the input cone mosaic.
    % In this visualization, the gray contours identify the mRGC RF centers
    hFig = figure(figNo); clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    hFig.Units = originalFigureUnits;
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
    if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(figureDir,'theMRGCMosaic.pdf'), hFig, 300)
    end
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
