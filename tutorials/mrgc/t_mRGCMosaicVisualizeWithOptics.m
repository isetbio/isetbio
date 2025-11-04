function t_mRGCMosaicVisualizeWithOptics(options)
% Visualize a prebaked mRGC mosaic and the optics used for its synthesis
%
% Description:
%   Demonstrates how to interactively load one of the available prebaked mRGC 
%   mosaics and visualize it, along with the PSF of used to synthesize it.
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

    t_mRGCMosaicVisualizeWithOptics()

    t_mRGCMosaicVisualizeWithOptics(...
        'croppedMosaicSizeDegs', [1 1]);

    t_mRGCMosaicVisualizeWithOptics(...
        'croppedMosaicSizeDegs', [2 1], ...
        'croppedMosaicEccentricityDegs', [-12 0]);

% ETTBSkip

    t_mRGCMosaicVisualizeWithOptics('promptUserForMosaic', true)

%}

arguments

    % Crop params
    options.croppedMosaicSizeDegs (1,:) double = [];
    options.croppedMosaicEccentricityDegs (1,:) double = [];

    % Whether to close previously open figures
    options.closePreviouslyOpenFigures (1,1) logical = true;

    % Prompt user for mosaic?
    options.promptUserForMosaic (1,1) logical = false;

    options.exportVisualizationPDF (1,1) logical = false;
    options.exportVisualizationPNG (1,1) logical = false;
end % arguments


% Set flags from key/value pairs

% Mosaic crop params 
croppedMosaicSizeDegs = options.croppedMosaicSizeDegs;
croppedMosaicEccentricityDegs = options.croppedMosaicEccentricityDegs;

exportVisualizationPDF = options.exportVisualizationPDF;
exportVisualizationPNG = options.exportVisualizationPNG;

% Close previously open figures
closePreviouslyOpenFigures = options.closePreviouslyOpenFigures;

if (closePreviouslyOpenFigures)
    % Close any stray figs
    close all;
end

% Obtain the directory where all the prebaked mRGCmosaics are stored
[theRGCmosaicFileNames, prebakedRGCmosaicDirectory] = mRGCMosaic.listPrebakedMosaics();

% Let the user select one prebaked mRGCmosaic
if (options.promptUserForMosaic)
    [mRGCMosaicFilename, prebakedMRGCMosaicDir] = ...
        uigetfile(sprintf('%s/*.mat', prebakedRGCmosaicDirectory), 'Select an RGB mosaic');
else
    prebakedMRGCMosaicDir = prebakedRGCmosaicDirectory;
    theFiles = dir(fullfile(prebakedMRGCMosaicDir,'*.mat'));
    if (isempty(theFiles))
        mRGCMosaicFilename = '';
    else
    mRGCMosaicFilename = theFiles(1).name;
end

% Load the user-selected prebaked mRGCmosaic
load(fullfile(prebakedMRGCMosaicDir,mRGCMosaicFilename), 'theMRGCMosaic');

% Optionally crop the mosaic
if (~isempty(croppedMosaicSizeDegs)) || (~isempty(croppedMosaicEccentricityDegs))
    theMRGCMosaic.cropToSizeAtEccentricity(...
            croppedMosaicSizeDegs, ...
            croppedMosaicEccentricityDegs);
end


% Generate the optics that were used to synthesize the mosaic 
visualizePSFonTopOfConeMosaic = ~true;
whichOptics = 'nativeOptics';
customRefractiveErrorDiopters = [];
[theOIatTheMosaicEccentricity, thePSFatTheMosaicEccentricity] = RGCMosaicAnalyzer.compute.opticsForResponses(...
        theMRGCMosaic, ...
        whichOptics, ...
        customRefractiveErrorDiopters, ...
        visualizePSFonTopOfConeMosaic);



% Get ready for publication-quality visualization
ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');


% Subdirectory for exporting the generated PDFs
exportVisualizationPDFdirectory = 'mosaicVisualizationPDFs';


% Visualize the full mosaic of RF centers using a representation
% like the representation used in visualizing 
% mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;


visualizedWidthDegs = theMRGCMosaic.sizeDegs(1);
visualizedHeightDegs = theMRGCMosaic.sizeDegs(2);
domainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 1];
domainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1];
domainVisualizationTicks = struct(...
    'x', theMRGCMosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 -0.5 0 0.5 1], ...
    'y', theMRGCMosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 -0.5 0 0.5 1]);

% Plot the mosaic of mRGC RF centers only
hFig = figure(1); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
ff.backgroundColor = [0 0 0];
ax = theAxes{1,1};

theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'identifyInputCones', false, ...
    'identifyPooledCones', false, ...
    'identifiedConeAperture', 'lightCollectingArea4sigma', ...
    'identifiedConeApertureThetaSamples', 16, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'centerSubregionContourSamples', 32, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
    'withFigureFormat', ff, ...
    'visualizationPDFfileName', sprintf('fullMRGCmosaicMinCenterConeWeight_%2.3f', minCenterConeWeight), ...
    'exportVisualizationPDF', true, ...
    'exportVisualizationPNG', true, ...
    'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);


% Visualize the full mosaic of RF centers, now visualizing the full
% extent of RF center pooling
hFig = figure(2); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
ax = theAxes{1,1};

% 1% min sensitivity for inclusion of divergent cone connections
minCenterConeWeight = mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections;

theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'identifyInputCones', false, ...
    'identifyPooledCones', false, ...
    'inputConesAlpha', 0.5, ...
    'identifiedConeAperture', 'lightCollectingArea4sigma', ...
    'identifiedConeApertureThetaSamples', 16, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'centerSubregionContourSamples', 32, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
    'withFigureFormat', ff, ...
    'visualizationPDFfileName', sprintf('fullMRGCmosaicMinCenterConeWeight_%2.3f', minCenterConeWeight), ...
    'exportVisualizationPDF', true, ...
    'exportVisualizationPNG', true, ...
    'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);



% Plot a smaller region of the mRGC mosaic with the PSF superimposed
narrowVisualizedWidthDegs = 0.2*min(theMRGCMosaic.sizeDegs(:));
narrowDomainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + [-1 1]*narrowVisualizedWidthDegs;
narrowDomainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + [-0.5 0.5]*narrowVisualizedWidthDegs;
narrowDomainVisualizationTicks = struct(...
    'x', -30:0.2:0, ...
    'y', -10:0.2:10);

% Generate a PSF visualization data struct (containing the vLambda-weighted PSF) for
% visualization purposes
PSFvisualizationOffset = theMRGCMosaic.eccentricityDegs - [mean(narrowDomainVisualizationLimits(1:2)) mean(narrowDomainVisualizationLimits(3:4))];
vLambdaWeightedPSF.data = RGCMosaicAnalyzer.compute.vLambdaWeightedPSF(thePSFatTheMosaicEccentricity);
vLambdaWeightedPSF.supportXdegs = thePSFatTheMosaicEccentricity.supportX/60 - PSFvisualizationOffset(1);
vLambdaWeightedPSF.supportYdegs = thePSFatTheMosaicEccentricity.supportY/60 - PSFvisualizationOffset(2);


hFig = figure(3); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
ax = theAxes{1,1};
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;
theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'identifyInputCones', true, ...
    'identifyPooledCones', true, ...
    'inputConesAlpha', 0.5, ...
    'identifiedConeAperture', 'lightCollectingArea4sigma', ...
    'identifiedConeApertureThetaSamples', 16, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'centerSubregionContourSamples', 32, ...
    'domainVisualizationLimits', narrowDomainVisualizationLimits, ...
    'domainVisualizationTicks', narrowDomainVisualizationTicks, ...
    'withSuperimposedPSF', vLambdaWeightedPSF, ...
    'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
    'withFigureFormat', ff, ...
    'visualizationPDFfileName', sprintf('zoomedInMRGCmosaicWithPSFminCenterConeWeight_%2.3f', minCenterConeWeight), ...
    'exportVisualizationPDF', exportVisualizationPDF, ...
    'exportVisualizationPNG', exportVisualizationPNG, ...
    'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);


end


