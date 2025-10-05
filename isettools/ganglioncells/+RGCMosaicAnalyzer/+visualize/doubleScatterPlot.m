%
% RGCMosaicAnalyzer.visualize.doubleScatterPlot(..)
%


function doubleScatterPlot(figNo, ...
	xData1, yData1, marker1, marker1Size, marker1Color, marker1FaceAlpha, marker1EdgeAlpha, populationContourInsteadOfPointCloud1, ...
	xData2, yData2, referenceYdata2, marker2, marker2Size, marker2Color, marker2FaceAlpha, marker2EdgeAlpha, populationContourInsteadOfPointCloud2, ...
	xLims, xTicks, yLims, yTicks, ...
	xScale, yScale, xAxisLabel, yAxisLabel, ...
	theLegends, thePDFFullFileName)

	% The density plot of acrom vs cone isolating BPIs
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    p1 = [];
    p2 = [];

    if (~isempty(xData1)) && (~isempty(yData1))
        if (~isempty(populationContourInsteadOfPointCloud1))
            xInterval = populationContourInsteadOfPointCloud1.xSupport(2)-populationContourInsteadOfPointCloud1.xSupport(1);
            xBinsNum = numel(populationContourInsteadOfPointCloud1.xSupport);
            yInterval = populationContourInsteadOfPointCloud1.ySupport(2)-populationContourInsteadOfPointCloud1.ySupport(1);
            yBinsNum = numel(populationContourInsteadOfPointCloud1.ySupport);

            theScatterPointDensityMap = zeros(yBinsNum, xBinsNum);

            for iPoint = 1:numel(xData1)
                xIndex = max([1 min([round(xData1(iPoint) / xInterval) xBinsNum])]);
                yIndex = max([1 min([round(yData1(iPoint) / yInterval) yBinsNum])]);
                theScatterPointDensityMap(yIndex, xIndex) = theScatterPointDensityMap(yIndex, xIndex) + 1;
            end

            % Normalize with respect to # of points in each xBin
            sums = sum(theScatterPointDensityMap,1);
            sums(sums == 0) = 1;
            theScatterPointDensityMap = bsxfun(@times, theScatterPointDensityMap, 1./sums);
            theScatterPointDensityMap = theScatterPointDensityMap / max(theScatterPointDensityMap(:));

            [X,Y] = meshgrid(...
                populationContourInsteadOfPointCloud1.xSupport, ...
                populationContourInsteadOfPointCloud1.ySupport+yInterval/2);

            for iP = 1:numel(populationContourInsteadOfPointCloud1.percentileLevels)
                lineStyle = 'none';
                theLevel = populationContourInsteadOfPointCloud1.percentileLevels(iP);
                contourf(ax, X,Y, theScatterPointDensityMap, theLevel*[1 1.01]/100, ...
                        'LineStyle', 'none', ...
                        'FaceColor', marker1Color*theLevel/100 + [1 1 1]*(1-theLevel/100), ...
                        'FaceAlpha', 0.8);
                hold(ax, 'on');
            end
            theLevel = populationContourInsteadOfPointCloud1.percentileLevels(1);
            [~,p1] = contour(ax, X,Y, theScatterPointDensityMap, theLevel*[1 1.01]/100, ...
                        'LineWidth', 1.0, 'LineStyle', '-', ...
                        'LineColor', marker1Color*0.5, ...
                        'EdgeAlpha', 0.4, ...
                        'FaceColor', 'none', ...
                        'FaceAlpha', 0.0);
            
        else
            p1 = scatter(ax, xData1, yData1, marker1Size^2, ...
            		'Marker', marker1, ...
            		'LineWidth', ff.lineWidth/2, ...
            		'MarkerFaceColor', marker1Color, ...
            		'MarkerEdgeColor', marker1Color*0.5, ...
            		'MarkerFaceAlpha', marker1FaceAlpha, ...
            		'MarkerEdgeAlpha', marker1EdgeAlpha);
        end
    end % (~isempty(xdata1)) && (~isempty(ydata1))

    if (~isempty(xData2)) && (~isempty(yData2))
        hold(ax, 'on');
        if (~isempty(populationContourInsteadOfPointCloud2))
            p2 = [];
        else
            p2 = scatter(ax, xData2, yData2, marker2Size^2, ...
        		'Marker', marker2, ...
        		'LineWidth', ff.lineWidth/2, ...
        		'MarkerFaceColor', marker2Color, ...
        		'MarkerEdgeColor', marker2Color*0.5, ...
        		'MarkerFaceAlpha', marker2FaceAlpha, ...
        		'MarkerEdgeAlpha', marker2EdgeAlpha);
        end

        if (~isempty(referenceYdata2))
            plot(ax, xData2, referenceYdata2, 'k--', 'LineWidth', 1.5);
        end

        set(ax, 'XScale', xScale, 'YScale', yScale);
        %axis (ax, 'square');

        set(ax, 'XLim', xLims, 'XTick', xTicks, 'YTick', yTicks);
        if (~isempty(yLims))
        	set(ax, 'YLim', yLims);
        else
        	yLims = get(ax, 'YLim');
        end
    end  % if (~isempty(xdata2)) && (~isempty(ydata2))


    
    if (~isempty(theLegends))
        if (~isempty(p1)) && (~isempty(p2))
    	   legendHandle = legend(ax, [p1 p2], theLegends, 'Orientation', 'Horizontal', 'NumColumns', 1, 'Location', 'NorthWest');
        else
            if (~isempty(p1))
                legendHandle = legend(ax, p1, theLegends{1}, 'Orientation', 'Horizontal', 'NumColumns', 1, 'Location', 'NorthWest');
            end
            if (~isempty(p2))
                legendHandle = legend(ax, p2, theLegends{2}, 'Orientation', 'Horizontal', 'NumColumns', 1, 'Location', 'NorthWest');
            end
        end
    	set(legendHandle, 'Color', [1 1 1], 'EdgeColor', 'none');
    end

	xlabel(ax, xAxisLabel);
	ylabel(ax, yAxisLabel);

	% Finalize figure using the Publication-Ready format
	ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, xLims, yLims);

    NicePlot.exportFigToPDF(thePDFFullFileName,hFig,  300);
end