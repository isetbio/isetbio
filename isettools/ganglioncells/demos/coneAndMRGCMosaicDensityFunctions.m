% Script to generate  mosaic density function plots which are use for some of
% the center-connectivity figures of our PLOS2024 paper
%

% Initialize session
close all; clear all;

theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();

hFig = RGCmodels.Watson.plot.figure1();
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'coneSpatialDensity.pdf'), hFig, 300);
    
hFig = RGCmodels.Watson.plot.figure9();
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'mRGCRFSpatialDensity.pdf'), hFig, 300);
    
hFig = RGCmodels.Watson.plot.figure14('inverseRatio', true);
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'ConeToMRGCRFDensityRatio.pdf'), hFig, 300);

