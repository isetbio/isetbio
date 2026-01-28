function t_mRGCMosaicValidationAgainstField2010(options)
% Visualize prebaked mRGCMosaics and the cone pooling maps of individual cells
%
% Description:
%    Validates against Field et al 2010:
%    (i)   how to load a pre-baked midget RGC mosaic, 
%    (ii)  visualize the mosaic of RF centers,
%    (iii) visualize the PSF of the optics used to build the mosaic
%    (iv)  visualize the RF center and surround cone pooling maps of
%          individual cells within the mosaic
%
%  This is set up with key/value pairs that demonstate how to select different
%  options. Different choices are illustrated in the examples
%  in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% History:
%    07/28/25  NPC  Wrote it.

% Examples:
%{

%
% NOTE: To run any RGC-related ISETBio code, such as this tutorial, users must follow
% the directions discribed in:
%    https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cell-(RGC)-mosaics
% under section "Configuring ISETBio to access RGC resources and run RGC simulations"
%

    t_mRGCMosaicValidationAgainstField2010(...
        'rgcMosaicName', 'PLOSpaperTemporal25DegsMosaic', ...
        'opticsSubjectName', 'PLOSpaperDefaultSubject', ...
        'syntheticMRGCpopulationEccMM', 6.75, ...
        'syntheticMRGCpopulationMarginMM', 0.3);


    t_mRGCMosaicValidationAgainstField2010(...
        'rgcMosaicName', 'PLOSpaperTemporal25DegsMosaic', ...
        'opticsSubjectName', 'PLOSpaperDefaultSubject', ...
        'syntheticMRGCpopulationEccMM', 6.25, ...
        'syntheticMRGCpopulationMarginMM', 0.3);

    t_mRGCMosaicValidationAgainstField2010(...
        'rgcMosaicName', 'PLOSpaperTemporal25DegsMosaic', ...
        'opticsSubjectName', 'PLOSpaperDefaultSubject', ...
        'syntheticMRGCpopulationEccMM', 6.00, ...
        'syntheticMRGCpopulationMarginMM', 0.3);

%}


arguments

    % ---- Mosaic specifiers for selecting a prebaked mRGC mosaic ---

    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'PLOSpaperNasal2DegsTinyMosaic';


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
            'PLOSpaperDefaultSubject' ...
            'PLOSpaperSecondSubject' ...
            'VSS2024TalkFirstSubject' ...
            'VSS2024TalkSecondSubject' ...
            'PLOSpaperStrehlRatio_0.87' ...
            'PLOSpaperStrehlRatio_0.72' ...
            'PLOSpaperStrehlRatio_0.59' ...
            'PLOSpaperStrehlRatio_0.60' ...
            'PLOSpaperStrehlRatio_0.27' ...
            'PLOSpaperStrehlRatio_0.23' ...
            'PLOSpaperStrehlRatio_0.21' ...
            'PLOSpaperStrehlRatio_0.19' ...
            'PLOSpaperStrehlRatio_0.09' ...
            } ...
            ) ...
        } ...
        = 'PLOSpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % See RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct
    % for all existing options
    options.targetVisualSTFdescriptor (1,:) char = 'default';


    % Ecc and margin (in mm) of synthetic mRGC population
    options.syntheticMRGCpopulationEccMM (1,1) double = 6.75
    options.syntheticMRGCpopulationMarginMM (1,1) double = 0.3

    % Whether to close previously open figures
    options.closePreviouslyOpenFigures (1,1) logical = true;
end


% Set flags from key/value pairs

% Mosaic specifiers for selecting a prebaked mRGC mosaic
rgcMosaicName = options.rgcMosaicName;
coneMosaicSpecies = options.coneMosaicSpecies;
opticsSubjectName = options.opticsSubjectName;
targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;

syntheticMRGCpopulationEccMM = options.syntheticMRGCpopulationEccMM;
syntheticMRGCpopulationMarginMM = options.syntheticMRGCpopulationMarginMM;
   

% Close previously open figures
closePreviouslyOpenFigures = options.closePreviouslyOpenFigures;
if (closePreviouslyOpenFigures)
    % Close any stray figs
    close all;
end


% Load the mRGCmosaic specified by the passed parameters:
%   coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor
% and generate the optics that were used to synthesize the mosaic
[theMRGCmosaic, ~, thePSFatTheMosaicEccentricity] = mRGCMosaic.loadPrebakedMosaic(...
        coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
        'computeTheMosaicOptics', true);


% Subdirectory for exporting the generated PDFs
exportVisualizationPDFdirectory = 'mosaicVisualizationPDFs';


% mRGC mosaic visualization limits and ticks (excluding the extent of the
% input cone mosaic)
visualizedWidthDegs = theMRGCmosaic.sizeDegs(1);
visualizedHeightDegs = theMRGCmosaic.sizeDegs(2);
domainVisualizationLimits(1:2) = theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 1];
domainVisualizationLimits(3:4) = theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1];
domainVisualizationTicks = struct(...
    'x', theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 -0.5 0 0.5 1], ...
    'y', theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 -0.5 0 0.5 1]);


% Visualize RF centers using the spatial extent of cones whose RF center pooling weights
% are >= mRGCMosaic.sensitivityAtPointOfOverlap;
% This representation is like the representation used in visualizing 
% mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;

% Get ready for publication-quality visualization
ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');
ff.backgroundColor = [0 0 0];

% Plot the mosaic of mRGC RF centers only
hFig = figure(1); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
ax = theAxes{1,1};

theMRGCmosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'identifyInputCones', false, ...
    'identifyPooledCones', false, ...
    'plottedRFoutlineFaceColor',  [0 1 0], ...
    'plottedRFoutlineFaceAlpha', 0.5, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'centerSubregionContourSamples', 32, ...
    'plottedRFoutlineFaceAlpha', 0.75, ...
    'plottedRFoutlineLineWidth', 1.0, ...
    'plottedRFoutlineFaceColor',  [0 1 0.4], ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
    'withFigureFormat', ff, ...
    'visualizationPDFfileName', sprintf('fullMRGCmosaicMinCenterConeWeight_%2.3f', minCenterConeWeight), ...
    'exportVisualizationPDF', true, ...
    'exportVisualizationPNG', true, ...
    'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);


% Visualize RF centers using the spatial extent of cones whose RF center pooling weights
% are > mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections
% (basically all cones connected to each RF center). This is useful for
% visualizing the degree of RF center overlap
minCenterConeWeight = mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections;

hFig = figure(2); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
ax = theAxes{1,1};

theMRGCmosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'identifyInputCones', false, ...
    'identifyPooledCones', false, ...
    'plottedRFoutlineFaceColor',  [0 1 0], ...
    'plottedRFoutlineFaceAlpha', 0.5, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'plottedRFoutlineFaceAlpha', 0.75, ...
    'plottedRFoutlineLineWidth', 1.0, ...
    'plottedRFoutlineFaceColor',  [0 1 0.4], ...
    'centerSubregionContourSamples', 32, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
    'withFigureFormat', ff, ...
    'visualizationPDFfileName', sprintf('fullMRGCmosaicMinCenterConeWeight_%2.3f', minCenterConeWeight), ...
    'exportVisualizationPDF', true, ...
    'exportVisualizationPNG', true, ...
    'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);





    % 10% which is the noise floor in the measurements of cone weights according to RF center overlap according to Greg Field
    % and which is what they used in Fig 4 of their 2010 paper.
    minCenterConeWeightForConePoolingMapVisualization = mRGCMosaic.minRFcenterConeWeightIncludedToMatchFigure4OfFieldEtAl2010;

    % Include surround cones whose pooling weights are >= 0.005 
    % (This is the threshold used in the Field et al 2010 paper -
    % Greg Field - personal communication)
    minSurroundConeWeightForConePoolingMapVisualization = 0.005;

    centerRadiiDegsX = zeros(1, theMRGCmosaic.rgcsNum);
    centerRadiiDegsY = zeros(1, theMRGCmosaic.rgcsNum);
    surroundRadiiDegs = zeros(1, theMRGCmosaic.rgcsNum);

    centerRadiiMicronsX = zeros(1, theMRGCmosaic.rgcsNum);
    centerRadiiMicronsY = zeros(1, theMRGCmosaic.rgcsNum);
    surroundRadiiMicrons = zeros(1, theMRGCmosaic.rgcsNum);

    eccDegs = zeros(1, theMRGCmosaic.rgcsNum);
    eccMicrons = zeros(1, theMRGCmosaic.rgcsNum);

    for iRGC = 1:theMRGCmosaic.rgcsNum
        [centerRadiiDegsX(iRGC), centerRadiiDegsY(iRGC), surroundRadiiDegs(iRGC), ...
         centerRadiiMicronsX(iRGC), centerRadiiMicronsY(iRGC), surroundRadiiMicrons(iRGC), ...
         eccDegs(iRGC), eccMicrons(iRGC)] = ...
            theMRGCmosaic.computeSpatialRFpoolingStats(iRGC, ...
                minCenterConeWeightForConePoolingMapVisualization, ...
                minSurroundConeWeightForConePoolingMapVisualization);
    end % iRGC

    % Only keep model cells with eccentricity are near the eccsyntheticMRGCpopulationEccMM
    idx = find(abs(eccMicrons(:)-syntheticMRGCpopulationEccMM*1000)<0.5*syntheticMRGCpopulationMarginMM*1000);

    centerRadiiDegsX = centerRadiiDegsX(idx);
    centerRadiiDegsY = centerRadiiDegsY(idx);
    surroundRadiiDegs = surroundRadiiDegs(idx);

    centerRadiiMicronsX = centerRadiiMicronsX(idx);
    centerRadiiMicronsY = centerRadiiMicronsY(idx);
    surroundRadiiMicrons = surroundRadiiMicrons(idx);

    eccDegs = eccDegs(idx);
    eccMicrons = eccMicrons(idx);




    minRadius = min([min(centerRadiiDegsX(:)) min(centerRadiiDegsY(:)) min(surroundRadiiDegs(:))]);
    maxRadius = max([max(centerRadiiDegsX(:)) max(centerRadiiDegsY(:)) max(surroundRadiiDegs(:))]);
    radiusBins = linspace(minRadius*2, maxRadius*2, 40);


    % Load the center/surround cone positions of the 3 cells in Field et al
    % (2010) Figure 4 E,F,G
    [cell1PoolingMap, cell2PoolingMap, cell3PoolingMap] = RGCmodels.Field2010.digitizedData.figure4PoolingMaps;


    % Fit an ellipse to the cones pooled by the RF center of cell1 to
    % compute the minor and major radii of the pool of cones pooled by the
    % RF center of cell1
    ellipseContourAngles = 0:10:350;
    maxNumberOfConesOutsideContour = 1;
    [~, ~, centerRadiusMicronsXField2010Cell1 , centerRadiusMicronsYField2010Cell1] = ...
        mRGCMosaic.subregionEllipseFromPooledConePositions(...
             cell1PoolingMap.centerConePositionsMicrons, [], ellipseContourAngles, [], maxNumberOfConesOutsideContour);

    % Compute the radius of the pool of cones pooled by the RF surround of cell1 
    surroundCenterMicrons = mean(cell1PoolingMap.centerConePositionsMicrons,1);
    surroundRadiusMicronsField2010Cell1 = max(sqrt(sum((bsxfun(@minus, cell1PoolingMap.surroundConePositionsMicrons, surroundCenterMicrons)).^2,2)));


    % Fit an ellipse to the cones pooled by the RF center of cell1 to
    % compute the minor and major radii of the pool of cones pooled by the
    % RF center of cell2
    [~, ~, centerRadiusMicronsXField2010Cell2 , centerRadiusMicronsYField2010Cell2] = ...
        mRGCMosaic.subregionEllipseFromPooledConePositions(...
             cell2PoolingMap.centerConePositionsMicrons, [], ellipseContourAngles, [], maxNumberOfConesOutsideContour);

    % Compute the radius of the pool of cones pooled by the RF surround of cell2
    surroundCenterMicrons = mean(cell2PoolingMap.centerConePositionsMicrons,1);
    surroundRadiusMicronsField2010Cell2 = max(sqrt(sum((bsxfun(@minus, cell2PoolingMap.surroundConePositionsMicrons, surroundCenterMicrons)).^2,2)));


    % Fit an ellipse to the cones pooled by the RF center of cell1 to
    % compute the minor and major radii of the pool of cones pooled by the
    % RF center of cell3
    [~, ~, centerRadiusMicronsXField2010Cell3 , centerRadiusMicronsYField2010Cell3] = ...
        mRGCMosaic.subregionEllipseFromPooledConePositions(...
             cell3PoolingMap.centerConePositionsMicrons, [], ellipseContourAngles, [], maxNumberOfConesOutsideContour);

    % Compute the radius of the pool of cones pooled b y the RF surround ofcell3
    surroundCenterMicrons = mean(cell3PoolingMap.centerConePositionsMicrons,1);
    surroundRadiusMicronsField2010Cell3 = max(sqrt(sum((bsxfun(@minus, cell3PoolingMap.surroundConePositionsMicrons, surroundCenterMicrons)).^2,2)));




    minRadius = min([min(centerRadiiMicronsX(:)) min(centerRadiiMicronsY(:)) min(surroundRadiiMicrons(:))]);
    maxRadius = max([max(centerRadiiMicronsX(:)) max(centerRadiiMicronsY(:)) max(surroundRadiiMicrons(:))]);
    radiusBins = 0:5:400;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(101); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    hold(ax, 'on');
    p1 = histogram(ax, centerRadiiMicronsX*2, radiusBins, 'FaceColor', [1 0.5 0.5]);
    p2 = histogram(ax, centerRadiiMicronsY*2, radiusBins, 'FaceColor', [1 0.7 0.15]);
    p3 = histogram(ax, surroundRadiiMicrons*2, radiusBins, 'FaceColor', [0.5 0.5 1.0]);

    yLims = get(ax, 'YLim');
    

    plot(ax, [centerRadiusMicronsXField2010Cell1*2 centerRadiusMicronsYField2010Cell1*2 surroundRadiusMicronsField2010Cell1*2], yLims(2)/2 * (1+0.1) *[1 1 1], 'k-', 'LineWidth', 1.5);
    plot(ax, [centerRadiusMicronsXField2010Cell1*2 centerRadiusMicronsYField2010Cell1*2 surroundRadiusMicronsField2010Cell1*2], yLims(2)/2 * (1) *[1 1 1], 'k-', 'LineWidth', 1.5);
    plot(ax, [centerRadiusMicronsXField2010Cell1*2 centerRadiusMicronsYField2010Cell1*2 surroundRadiusMicronsField2010Cell1*2], yLims(2)/2 * (1-0.1) *[1 1 1], 'k-', 'LineWidth', 1.5);

    scatter(ax, centerRadiusMicronsXField2010Cell1*2, yLims(2)/2 * (1+0.1), 144, 'o', ...
        'MarkerFaceColor', [1 0.5 0.5], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    scatter(ax, centerRadiusMicronsXField2010Cell2*2, yLims(2)/2, 144, 'o', ...
        'MarkerFaceColor', [1 0.5 0.5], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    scatter(ax, centerRadiusMicronsXField2010Cell3*2, yLims(2)/2* (1-0.1), 144, 'o', ...
        'MarkerFaceColor', [1 0.5 0.5], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    xlabel(ax, 'subregion diameter (microns)');

    
    scatter(ax, centerRadiusMicronsYField2010Cell1*2, yLims(2)/2* (1+0.1), 144, 's', ...
        'MarkerFaceColor', [1 0.7 0.15], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    scatter(ax, centerRadiusMicronsYField2010Cell2*2, yLims(2)/2, 144, 's', ...
        'MarkerFaceColor', [1 0.7 0.15], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
     scatter(ax, centerRadiusMicronsYField2010Cell3*2, yLims(2)/2* (1-0.1), 144, 's', ...
        'MarkerFaceColor', [1 0.7 0.15], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);


    
    scatter(ax, surroundRadiusMicronsField2010Cell1*2, yLims(2)/2* (1+0.1), 200, 'h', ...
        'MarkerFaceColor', [0.5 0.5 1.0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    scatter(ax, surroundRadiusMicronsField2010Cell2*2, yLims(2)/2, 200, 'h', ...
        'MarkerFaceColor', [0.5 0.5 1.0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
    scatter(ax, surroundRadiusMicronsField2010Cell3*2, yLims(2)/2* (1-0.1), 200, 'h',...
        'MarkerFaceColor', [0.5 0.5 1.0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'LineWidth', 1.0);
   
    


    ylabel(ax, '# of cells');
    set(ax, 'XLim', [0 400], 'XTick', [0 50 100 150 200 250 300 350 400], 'XTickLabel', {'0', '', '100', '', '200', '', '300', '', '400'});

    legend(ax, [p1 p2 p3], {'center (major)', 'center (minor)', 'surround'});

    pdfRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theDir = fullfile(pdfRootDir, exportVisualizationPDFdirectory);

    % The PDF filename
    thePDFFileName = sprintf('%s/ConePoolingStats%2.2fMM.pdf',theDir, syntheticMRGCpopulationEccMM);
   
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);
    drawnow;
    NicePlot.exportFigToPDF(thePDFFileName, hFig, 300);



end



