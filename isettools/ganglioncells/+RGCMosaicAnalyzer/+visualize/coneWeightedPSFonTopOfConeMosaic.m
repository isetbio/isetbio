function coneWeightedPSFonTopOfConeMosaic(figNo, theInputConeMosaic, theConeWeightedPSF, thePDFfileName, theTitle)
% CONEWEIGHTEDPSFONTOPOFCONEMOSAIC Visualizes the cone-weighted Point Spread Function (PSF)
% on top of a cone mosaic.
%
% (Written by CoPilot)
%
% Syntax:
%   coneWeightedPSFonTopOfConeMosaic(figNo, theInputConeMosaic, theConeWeightedPSF,
%   thePDFfileName, theTitle)
%
% Description:
%   This function creates a visualization of the cone-weighted Point Spread Function (PSF)
%   overlaid on a specified cone mosaic. It generates a figure with axes for displaying the
%   PSF and exports the visualization to a PDF file.
%
% Inputs:
%   figNo                    - Figure number for the visualization.
%   theInputConeMosaic       - The cone mosaic object to visualize the PSF on.
%   theConeWeightedPSF       - The cone-weighted PSF data to visualize.
%   thePDFfileName           - Name of the PDF file to export the visualization.
%   theTitle                 - Title for the visualization plot.
%
% Example:
%   coneWeightedPSFonTopOfConeMosaic(1, coneMosaic, coneWeightedPSF, 'output.pdf', 'Cone PSF Visualization');
%
% See also: VISUALIZEPSFATWAVELENGTHONTOPOFCONEMOSAIC, NICEPLOT

ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
hFig = figure(figNo); clf;
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

ax = theAxes{1,1};
visualizedWavelengthIndex = 1;
maxPSF = 1;
infoString = '';
xo = theInputConeMosaic.eccentricityDegs(1);
yo = theInputConeMosaic.eccentricityDegs(2);
RGCMosaicConstructor.helper.optics.visualizePSFAtWavelengthOnTopOfConeMosaic(...
    hFig, ax, theInputConeMosaic, theConeWeightedPSF, ...
    visualizedWavelengthIndex, theConeWeightedPSF, ...
    theTitle, ...
    'XLimsArcMin', [-12.5 12.5], ...
    'YLimsArcMin', [-12.5 12.5], ...
    'XTicksArcMin', xo*60 + (-12:6:12), ...
    'YTicksArcMin', yo*60 + (-12:6:12), ...
    'withInfoString', infoString);

ff.box = 'on';
ff.tickDir = 'in';
xlabel(ax, 'eccentricity, x (degs)');
ylabel(ax, 'eccentricity, y (degs)');
% Finalize figure using the Publication-Ready format
PublicationReadyPlotLib.applyFormat(ax,ff);
%PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

% OLD WAY
%theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();

%if (isdir(theRawFiguresDir))
%    pdfExportSubDir = 'validation';

%    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFfileName);
%    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
%else
%    fprintf('Raw figures directory not found: %s.\nWill not export PDF\n', theRawFiguresDir);
%end


p = getpref('isetbio');
pdfExportSubDir = fullfile(p.rgcResources.figurePDFsDir);
theVisualizationPDFfilename = fullfile('modelComponentVisualizations', thePDFfileName);

% Generate the path if we need to
RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
    pdfExportSubDir, theVisualizationPDFfilename, ...
    'generateMissingSubDirs', true);

thePDFfileName = fullfile(pdfExportSubDir, theVisualizationPDFfilename);
NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);

end