function achromaticAndConeIsolatingSTFs(figNo, theRGCindex, rfCenterConeDominance, ...
			theSpatialFrequencySupport, theAchromaticSTF, ...
			theLconeIsolatingSTF, theMconeIsolatingSTF, ...
			theAcromaticSTFBPI, theLconeIsolatingSTFBPI, theMconeIsolatingSTFBPI,...
			maxSTF, pdfFileName, varargin)

	% Parse input
    p = inputParser;
    % Optional params
    p.addParameter('spatialFrequencyLims', [0.05 150], @(x)(isnumeric(x)&&(numel(x)==2)));
    p.addParameter('squareAxes', true, @islogical);
    p.addParameter('lineWidthBoost', 1.0, @isscalar);
    p.addParameter('markerSizeBoost', 0.0, @isscalar);
    p.parse(varargin{:});
	XLims = p.Results.spatialFrequencyLims;
	squareAxes = p.Results.squareAxes;
	lineWidthBoost = p.Results.lineWidthBoost;
	markerSizeBoost = p.Results.markerSizeBoost;

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = 'validation';

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    LcenterColor = RGCMosaicConstructor.constants.LcenterColor;
	McenterColor = RGCMosaicConstructor.constants.McenterColor;
	achromaticColor = RGCMosaicConstructor.constants.achromaticColor;

	plot(ax, theSpatialFrequencySupport, theAchromaticSTF, '-', ...
		'LineWidth', 6, 'Color', achromaticColor*0);
	hold(ax, 'on');
	plot(ax, theSpatialFrequencySupport, theAchromaticSTF, '-', ...
		'LineWidth', 3, 'Color', achromaticColor);

	plot(ax, theSpatialFrequencySupport, theAchromaticSTF, 'o', 'MarkerSize', ff.markerSize+markerSizeBoost, ...
		'MarkerFaceColor', achromaticColor, 'MarkerEdgeColor', achromaticColor*0, ...
		'LineWidth', ff.lineWidth*lineWidthBoost, 'Color', achromaticColor*0);
	
	if (rfCenterConeDominance == cMosaic.LCONE_ID)
		plot(ax, theSpatialFrequencySupport, theLconeIsolatingSTF, '-', ...
			'LineWidth', 6, 'Color', LcenterColor*0);
		plot(ax, theSpatialFrequencySupport, theLconeIsolatingSTF, '-', ...
			'LineWidth', 3, 'Color', LcenterColor);

		plot(ax, theSpatialFrequencySupport, theLconeIsolatingSTF, 'o', 'MarkerSize', ff.markerSize+markerSizeBoost, ...
			'MarkerFaceColor', LcenterColor, 'MarkerEdgeColor', LcenterColor*0, ...
			'LineWidth', ff.lineWidth*lineWidthBoost, 'Color', LcenterColor*0);
	end

	if (rfCenterConeDominance == cMosaic.MCONE_ID)
		plot(ax, theSpatialFrequencySupport, theMconeIsolatingSTF, '-', ...
			'LineWidth', 6, 'Color', McenterColor*0);
		plot(ax, theSpatialFrequencySupport, theMconeIsolatingSTF, '-', ...
			'LineWidth', 3, 'Color', McenterColor);

		plot(ax, theSpatialFrequencySupport, theMconeIsolatingSTF, 'o', 'MarkerSize', ff.markerSize+markerSizeBoost, ...
			'MarkerFaceColor', McenterColor, 'MarkerEdgeColor', McenterColor*0, ...
			'LineWidth', ff.lineWidth*lineWidthBoost, 'Color', LcenterColor*0);
	end

	if (squareAxes)
    	axis (ax, 'square');
    end

    XTicks = [0.03 0.1 0.3 1 3 10 30 100 150];
    XTickLabels =  {'.03', '.1', '.3', '1', '3', '10', '30', '100', ''};

    YLims = [0 ceil(maxSTF*10)/10];
    set(ax, 'XTick', XTicks, ...
    	'XTickLabel', XTickLabels, ...
    	'YTick', 0:0.1:1);
    set(ax, 'XScale', 'log', 'XLim', XLims, 'YLim', YLims);

    xlabel(ax,'spatial frequency (c/deg)');
    ylabel(ax,'response modulation');
    title(ax, sprintf('RGC #%d, BPIs: %2.2f %2.2f %2.2f', theRGCindex, theAcromaticSTFBPI, theLconeIsolatingSTFBPI, theMconeIsolatingSTFBPI));

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);

end
