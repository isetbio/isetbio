%
% RGCMosaicConstructor.visualize.spectralUniformitySpectralCompactnessHistograms(theMRGCMosaic)
%
function spectralUniformitySpectralCompactnessHistograms(theMRGCMosaic)

	% Compute thePatchSpectralUniformity, thePatchCentroidOverlap, thePatchCentroidSigma
     [theSpatialCompactnessCosts, theSpectralUniformityCosts, ~, ~, ...
         theCentroidOverlapCosts, theSpatialVarianceCosts] = theMRGCMosaic.rfCenterSpatioChromaticCosts();

     theRGCConeMixtures = 0.5 + (0.5-theSpectralUniformityCosts/2);
     coneMixtureTickLabels = {'', '50:50', '', '75:25', '', '100:0'};
	theConeMixtureBins = [0.5 0.6 0.7 0.8 0.9 1.0];

     theRGCCentroidNormalizedDistances = 1./(2*theCentroidOverlapCosts);
     theRGCCentroidNormalizedDistanceBins = 0:0.125:(1.5);
     theRGCCentroidNormalizedDistanceLabels = {'0', '',  '.25', '', '.5', '', '.75', '', '1',  '',  '1.25', '', '>=1.5'};

     % The means
     meanDistance = mean(theRGCCentroidNormalizedDistances(:));
     meanConeMixture = mean(theRGCConeMixtures(:));

	faceColor = [1 0.75 0.5];
	edgeColor = [1 0.5 0];
	faceAlpha = 0.75;
	lineWidth = 1.5;

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir()
	pdfExportSubDir = 'centerConnected';
	
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(10); clf;
     theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
     ax = theAxes{1,1};
     [theCounts, theBarPlotHandle] = RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(...
          ax, theRGCCentroidNormalizedDistances, theRGCCentroidNormalizedDistanceBins, ...
    		faceColor, edgeColor, faceAlpha, lineWidth, ...
    		'binWidth', []);
     hold(ax, 'on');
     

     plot(ax, meanDistance*[1 1], [0 1], 'k--', 'LineWidth', 1.5);
     XLims = [0.25-0.125 1.5+0.125/2];
     YLims = [0 1]; YTick = 0:0.2:1.0;
	set(ax, 'XLim', XLims, 'XTick', theRGCCentroidNormalizedDistanceBins, 'XTickLabel', theRGCCentroidNormalizedDistanceLabels, ...
		'YLim', YLims, 'YTick', YTick);
     title(ax, sprintf('mean: %2.3f', meanDistance));
     xlabel(ax,'RGC centroid normalized distance');
     ylabel(ax, 'frequency')

     % Finalize figure using the Publication-Ready format
     ff.box = 'off';
     PublicationReadyPlotLib.applyFormat(ax,ff);

     % Export figure
     theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
     thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'spatialCompactnessHistogram.pdf');
     NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);




	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(30); clf;
     theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
     ax = theAxes{1,1};

     [theCounts, theBarPlotHandle] = RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(...
          ax, theRGCConeMixtures, theConeMixtureBins, ...
    		faceColor, edgeColor, faceAlpha, lineWidth, ...
    		'binWidth', []);

     hold(ax, 'on');
     plot(ax, meanConeMixture*[1 1], [0 1], 'k--', 'LineWidth', 1.5);
     XLims = [0.5 1.101];
     YLims = [0 1]; YTick = 0:0.2:1.0;
	set(ax, 'XLim', XLims, 'XTick', theConeMixtureBins, 'XTickLabel', coneMixtureTickLabels, ...
		'YLim', YLims, 'YTick', YTick);
     title(ax, sprintf('mean: %2.3f', meanConeMixture));
     xlabel(ax,'dominant:non-dominant cone mixture');
     ylabel(ax, 'frequency')

     % Finalize figure using the Publication-Ready format
     ff.box = 'off';
     PublicationReadyPlotLib.applyFormat(ax,ff);

     % Export figure
     theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
     thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'spectralUniformityHistogram.pdf');
     NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end
