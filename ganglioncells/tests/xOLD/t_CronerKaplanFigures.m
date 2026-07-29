% Script to generate some Croner & Kaplan data figures fof our PLOS2024 paper
%
% t_CronerKaplanFigures
%

% Initialize session
close all; clear all;

theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();

figNo = 10;
[hFig, hFig2, hFig3] = RGCmodels.CronerKaplan.plot.typicalSTF(figNo);
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'macaqueRGCtypicalSTF.pdf'), hFig, 300);
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'macaqueRGCtypicalRsRc.pdf'), hFig2, 300);
NicePlot.exportFigToPDF(fullfile(theRawFiguresDir, 'macaqueRGCtypicalIntSCsens.pdf'), hFig3, 300);
