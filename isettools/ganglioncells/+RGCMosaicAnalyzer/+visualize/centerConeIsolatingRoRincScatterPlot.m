function centerConeIsolatingRoRincScatterPlot(figNo, ...
	achromaticRatios, centerConeIsolatingRatios, centerConeDominances, ...
	pdfFileName, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('markerSize', 12, @isscalar);
    p.addParameter('markerFaceAlpha', 0.6, @isscalar);
    p.addParameter('markerLineWidth', 1.0, @isscalar);
    p.addParameter('markerEdgeAlpha', 0.5, @isscalar);
    p.addParameter('showLegends', false, @islogical);
	p.parse(varargin{:});

	markerSize = p.Results.markerSize;
	markerFaceAlpha = p.Results.markerFaceAlpha;
	markerEdgeAlpha = p.Results.markerEdgeAlpha;
	markerLineWidth = p.Results.markerLineWidth;
	showLegends = p.Results.showLegends;

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo+1); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    hold(ax, 'on');

    idx = find(centerConeDominances == cMosaic.LCONE_ID);
    centerColor = RGCMosaicConstructor.constants.LcenterColor;
    scatter(ax, achromaticRatios(idx), centerConeIsolatingRatios(idx), markerSize^2, ...
    	'MarkerFaceColor', centerColor, 'MarkerEdgeColor', [0 0 0], ...
    	'MarkerFaceAlpha', markerFaceAlpha, ...
    	'MarkerEdgeAlpha', markerEdgeAlpha, ...
    	'LineWidth', markerLineWidth);

    idx = find(centerConeDominances == cMosaic.MCONE_ID);
    centerColor = RGCMosaicConstructor.constants.McenterColor;
    scatter(ax, achromaticRatios(idx), centerConeIsolatingRatios(idx), markerSize^2, ...
    	'MarkerFaceColor', centerColor, 'MarkerEdgeColor', [0 0 0], ...
    	'MarkerFaceAlpha', markerFaceAlpha, ...
    	'MarkerEdgeAlpha', markerEdgeAlpha, ...
    	'LineWidth', markerLineWidth);

    axis (ax, 'square');
    XLims = [0 1.0]; YLims = [0 1.0];
	grid(ax, 'on');
	set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
		    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
    xlabel(ax,'achromatic RF map');
    ylabel(ax,'center cone-isolating RF map');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);


    if (~isempty(pdfFileName))
    	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    	pdfExportSubDir = 'validation';
    	thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end
end
