function [purityBins, theLcenterCounts, theMcenterCounts] = surroundPurityDistribution(surroundConePurities, radialEccentricities, ...
	LcenterCellIndices, McenterCellIndices, visualizedRadialEccentricityRange)


	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir()
	pdfExportSubDir = 'validation';
	
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(10); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    purityBins = 0:0.02:1.0;
    faceAlpha = 0.6;
    lineWidth = 1.5;

    theLcenterCounts = RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, surroundConePurities(LcenterCellIndices), purityBins, ...
    	 RGCMosaicConstructor.constants.LcenterColor, RGCMosaicConstructor.constants.LcenterColor*0.5, faceAlpha, lineWidth, ...
    	 'relativeHeight', 1.0);
    hold(ax, 'on');

    theMcenterCounts = RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, surroundConePurities(McenterCellIndices), purityBins, ...
    	 RGCMosaicConstructor.constants.McenterColor, RGCMosaicConstructor.constants.McenterColor*0.5, faceAlpha, lineWidth, ...
    	 'relativeHeight', numel(McenterCellIndices)/numel(LcenterCellIndices));

	XLims = [0 1];
	YLims = get(ax, 'YLim');
	set(ax, 'XLim', XLims, 'XTick', [0 1/6 1/3 1/2 2/3 5/6 1], 'XTickLabel', {'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'}, ...
		'YLim', YLims, 'YTickLabel', {});
	legend(ax, {'L-center', 'M-center'});
    xlabel(ax,'surround cone purity');
    ylabel(ax, 'frequency')
    %title(ax,'Sp(L-center)=wS(M)/(wS(L)+wS(M), Sp(M-center)=S(L)/(wS(L)+wS(M)');
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'population_surround_cone_purity.pdf');
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);


    hFig = figure(11); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    scatter(ax,radialEccentricities(LcenterCellIndices), surroundConePurities(LcenterCellIndices), 100, ...
    	'MarkerFaceColor', RGCMosaicConstructor.constants.LcenterColor, 'MarkerEdgeColor', RGCMosaicConstructor.constants.LcenterColor*0.5, 'MarkerFaceAlpha', 0.4);
    hold(ax, 'on');
    scatter(ax,radialEccentricities(McenterCellIndices), surroundConePurities(McenterCellIndices), 100, ...
    	'MarkerFaceColor', RGCMosaicConstructor.constants.McenterColor, 'MarkerEdgeColor', RGCMosaicConstructor.constants.McenterColor*0.5, 'MarkerFaceAlpha', 0.4);
    XLims = visualizedRadialEccentricityRange;
    YLims = [0 1];
    set(ax, 'YLim', YLims, 'YTick', [0 1/6 1/3 1/2 2/3 5/6 1], 'YTickLabel', {'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'}, ...
    	'XLim', XLims, 'XTick', visualizedRadialEccentricityRange(1):0.5:visualizedRadialEccentricityRange(2));
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'surround cone purity');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

   
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'surround_cone_purity_vs_eccentricity.png');
    NicePlot.exportFigToPNG(thePDFfileName,hFig,  300);
end

