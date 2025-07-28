function LeeShapleyBPIdata()
%
% RGCMosaicAnalyzer.visualize.LeeShapleyBPIdata()
%

    dataOut = RGCMosaicConstructor.publicData.LeeShapley.midgetBandPassIndices();
    centerColor = RGCMosaicConstructor.constants.LcenterColor;

    % Only midget RGC cells
    LeeShapleyLcenterBPIdata = dataOut('L-centerMidget');
    LeeShapleyMcenterBPIdata = dataOut('M-centerMidget');

    % Only parvo LGN cells
    LeeShapleyLcenterBPIdata = dataOut('L-centerParvo');
    LeeShapleyMcenterBPIdata = dataOut('M-centerParvo');


    % Combined midget and parvo cells
    LeeShapleyLcenterBPIdata = dataOut('L-center');
    LeeShapleyMcenterBPIdata = dataOut('M-center');

    showLcenterData = true;
    showMcenterData = true;

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = 'validation';

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(222); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    if (showLcenterData)
        scatter(ax, LeeShapleyLcenterBPIdata.x, LeeShapleyLcenterBPIdata.y, ff.markerSize^2, ...
    		'Marker', 's', ...
    		'MarkerFaceColor', RGCMosaicConstructor.constants.LcenterColor, ...
    		'MarkerEdgeColor', [0 0 0], ...
    		'MarkerFaceAlpha', 0.5);
    end

    hold(ax, 'on');
    if (showMcenterData)
        scatter(ax, LeeShapleyMcenterBPIdata.x, LeeShapleyMcenterBPIdata.y, ff.markerSize^2, ...
    		'Marker', 's', ...
    		'MarkerFaceColor', RGCMosaicConstructor.constants.McenterColor, ...
    		'MarkerEdgeColor', [0 0 0], ...
    		'MarkerFaceAlpha', 0.5);
    end

    %plot(ax, [0 1], [0.5 1], 'k-', 'LineWidth', 2.0);

    theLegends = {};
    if (showLcenterData)
        theLegends{numel(theLegends)+1} = 'L-center'
    end
    if (showMcenterData)
        theLegends{numel(theLegends)+1} = 'M-center'
    end

    if (~isempty(theLegends))
        legend(ax, theLegends, 'Location', 'SouthWest');
    end

    axis (ax, 'square');
    XLims = [0 1.0]; YLims = [0 1.0];
	grid(ax, 'on');
	set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
		    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
    xlabel(ax, 'achromatic gratings');
    ylabel(ax, 'cone-isolating gratings');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    pdfFileName = 'LeeShapleyRawData.pdf'
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end
