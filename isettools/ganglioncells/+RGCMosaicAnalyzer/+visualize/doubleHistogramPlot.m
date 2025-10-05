%
% RGCMosaicAnalyzer.visualize.doubleHistogramPlot(..)
%

function doubleHistogramPlot(figNo, ...
	data1, bins1, color1, binWidth1, ...
	data2, bins2, color2, binWidth2, ...
	xLims, xTicks, yLims, xAxisLabel, yAxisLabel, ...
	theLegends, thePDFFullFileName)

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, data1, bins1, ...
    	 color1, color1*0.5, 0.7, ff.lineWidth, ...
    	 'relativeHeight', 1.0, ...
    	 'binWidth', binWidth1);
    hold(ax, 'on');
    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, data2, bins2, ...
    	 color2, color2*0.5, 0.7, ff.lineWidth, ...
    	 'relativeHeight', 1.0, ...
    	 'binWidth', binWidth2);

    %axis(ax, 'square');
	set(ax, 'XLim', xLims, 'XTick', xTicks, 'YTick', 0:0.1:1); 

	
	if (~isempty(yLims))
    	set(ax, 'YLim', yLims);
    else
    	yLims = get(ax, 'YLim');
    end

    if (~isempty(theLegends))
    	legendHandle = legend(ax, theLegends, 'Orientation', 'Horizontal', 'NumColumns', 1, 'Location', 'NorthEast');
    	set(legendHandle, 'Color', [1 1 1], 'EdgeColor', 'none');
	end

    xlabel(ax, xAxisLabel);
    ylabel(ax, yAxisLabel)

    % Finalize figure using the Publication-Ready format
    ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, xLims, yLims);

    NicePlot.exportFigToPDF(thePDFFullFileName,hFig,  300);
end