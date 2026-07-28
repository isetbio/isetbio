function scatterPrcTilesPlot(figNo, ...
		xData1, yPrctTileData1, marker1, marker1Size, marker1Color, marker1FaceAlpha, marker1EdgeAlpha, ...
	    xData2, yData2, referenceYdata2, marker2, marker2Size, marker2Color, marker2FaceAlpha, marker2EdgeAlpha, ...
	    xLims, xTicks, yLims, yTicks, ...
	    xScale, yScale, xAxisLabel, yAxisLabel, theLegends, ...
        exportVisualizationPDF, thePDFFullFileName)

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    p1 = [];
    p2 = [];

  
    if (~isempty(xData2)) && (~isempty(yData2))
        
        p2 = scatter(ax, xData2, yData2, marker2Size^2, ...
        		'Marker', marker2, ...
        		'LineWidth', ff.lineWidth/2, ...
        		'MarkerFaceColor', marker2Color, ...
        		'MarkerEdgeColor', marker2Color*0.5, ...
        		'MarkerFaceAlpha', marker2FaceAlpha, ...
        		'MarkerEdgeAlpha', marker2EdgeAlpha);
    
        hold(ax, 'on');
        if (~isempty(referenceYdata2))
            plot(ax, xData2, referenceYdata2, 'k-', 'LineWidth', 3);
            plot(ax, xData2, referenceYdata2, 'w-', 'LineWidth', 1.5);
        end
    
        set(ax, 'XScale', xScale, 'YScale', yScale);
    
        set(ax, 'XLim', xLims, 'XTick', xTicks, 'YTick', yTicks);
        if (~isempty(yLims))
        	set(ax, 'YLim', yLims);
        else
        	yLims = get(ax, 'YLim');
        end
    end

    if (~isempty(xData1)) && (~isempty(yPrctTileData1))

        for i = 1:numel(xData1)
            xData(i) = mean(xData1{i});
        end

        hold(ax, 'on');

        %plot(ax, xData, yPrctTileData1(:,2), 'k-', ...
        %    		'LineWidth', ff.lineWidth);
        
        lowestPrcTileIndex = 1;
        maxPrcTileIndex = size(yPrctTileData1,2);
        yDataUpper = yPrctTileData1(:,lowestPrcTileIndex);
        yDataLower = yPrctTileData1(:,maxPrcTileIndex);

        v = [];
        for i = 1:numel(xData)
            v = cat(1,v,  [xData(i) yDataUpper(i)]);
        end
        for i = numel(xData):-1:1
            v = cat(1,v, [xData(i) yDataLower(i)]);
        end

        f = 1:size(v,1);

        p1 = patch('Faces',f,'Vertices',v,...
            'EdgeColor',marker1Color*0.5,'FaceColor',marker1Color,'FaceAlpha', marker1FaceAlpha, 'EdgeAlpha', 0.0, 'LineWidth',ff.lineWidth/2);


        for kk = 1:2
            yDataUpper = yPrctTileData1(:,lowestPrcTileIndex+kk);
            yDataLower = yPrctTileData1(:,maxPrcTileIndex-kk);

            v = [];
            for i = 1:numel(xData)
                v = cat(1,v,  [xData(i) yDataUpper(i)]);
            end
            for i = numel(xData):-1:1
                v = cat(1,v, [xData(i) yDataLower(i)]);
            end

            f = 1:size(v,1);

            patch('Faces',f,'Vertices',v,...
                'EdgeColor',marker1Color*0.5,'FaceColor',marker1Color,'FaceAlpha', marker1FaceAlpha, 'EdgeAlpha', 0.0, 'LineWidth',ff.lineWidth/2);

        end

        if (1==2)
        p1 = scatter(ax, xData, yPrctTileData1(:,2), marker1Size^2, ...
            		'Marker', marker1, ...
            		'LineWidth', ff.lineWidth/2, ...
            		'MarkerFaceColor', marker1Color, ...
            		'MarkerEdgeColor', marker1Color*0.5, ...
            		'MarkerFaceAlpha', marker1FaceAlpha, ...
            		'MarkerEdgeAlpha', marker1EdgeAlpha);
        end
        

    end

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

    if (exportVisualizationPDF)
        NicePlot.exportFigToPDF(thePDFFullFileName,hFig,  300, 'beVerbose');
    end
end