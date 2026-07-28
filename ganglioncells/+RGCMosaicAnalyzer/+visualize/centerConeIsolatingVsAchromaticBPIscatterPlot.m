function centerConeIsolatingVsAchromaticBPIscatterPlot(figNo, ...
	achromaticBPIs, coneIsolatingBPIs, surroundConePurities, centerConeDominance, ...
	pdfFileName, pdfFileName2, pdfFileName3, pdfFileName4, pdfFileName5, varargin)
	
	% Parse input
    p = inputParser;
    p.addParameter('markerSize', 12, @isscalar);
    p.addParameter('markerFaceAlpha', 0.3, @isscalar);
    p.addParameter('markerLineWidth', 1.0, @isscalar);
    p.addParameter('markerEdgeAlpha', 1.0, @isscalar);
    p.addParameter('onlyDepictLeeShapleyData', false, @islogical);
    p.addParameter('superimposeLeeShapleyData', false, @islogical);
    p.addParameter('showLegends', false, @islogical);
    p.addParameter('superimposeLeeShapleyDataMarkerSize', 14, @isscalar);
    p.addParameter('maxVisualizedBPIDensity', 1.0, @isscalar);
    p.addParameter('bpiInterval', 0.025, @isscalar);
    p.addParameter('minNumberOfCellsPerBPIintervalForPercentilePlot', 2, @isscalar);
    p.addParameter('thePercentileValues', [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95], @isnumeric);

    p.parse(varargin{:});
	markerSize = p.Results.markerSize;
	markerFaceAlpha = p.Results.markerFaceAlpha;
	markerEdgeAlpha = p.Results.markerEdgeAlpha;
	markerLineWidth = p.Results.markerLineWidth;
	onlyDepictLeeShapleyData = p.Results.onlyDepictLeeShapleyData;
	superimposeLeeShapleyData = p.Results.superimposeLeeShapleyData;
	superimposeLeeShapleyDataMarkerSize = p.Results.superimposeLeeShapleyDataMarkerSize;
	maxVisualizedBPIDensity = p.Results.maxVisualizedBPIDensity;
	bpiInterval = p.Results.bpiInterval;
	minNumberOfCellsPerBPIintervalForPercentilePlot = p.Results.minNumberOfCellsPerBPIintervalForPercentilePlot;
	thePercentileValues = p.Results.thePercentileValues;
	showLegends = p.Results.showLegends;

	bpiBins = 0:bpiInterval:(1.0-bpiInterval);
	maxBPIindex = numel(bpiBins);
	bpiMap = zeros(maxBPIindex,maxBPIindex);

	histData.achromaticBPIs = bpiBins;
	histData.coneIsolatingBPIs = cell(1, numel(histData.achromaticBPIs));
	for iBin = 1:numel(histData.achromaticBPIs)
		histData.coneIsolatingBPIs{iBin} = [];
	end
	

	for iRGC = 1:numel(achromaticBPIs)
		achromaticBPIindex = max([1 min([round(achromaticBPIs(iRGC)/bpiInterval) maxBPIindex])]);
		coneIsolatingBPIindex = max([1 min([round(coneIsolatingBPIs(iRGC)/bpiInterval) maxBPIindex])]);
		d = histData.coneIsolatingBPIs{achromaticBPIindex};
		d(numel(d)+1) = coneIsolatingBPIs(iRGC);
		histData.coneIsolatingBPIs{achromaticBPIindex} = d;
		bpiMap(coneIsolatingBPIindex, achromaticBPIindex) = bpiMap(coneIsolatingBPIindex, achromaticBPIindex) + 1;
	end

	
	for iBin = 1:numel(histData.achromaticBPIs)
		d = histData.coneIsolatingBPIs{iBin};
		if (numel(d) > minNumberOfCellsPerBPIintervalForPercentilePlot)
			histData.prctileConeIsolatingBPIindex(iBin,:) = prctile(d(:), thePercentileValues);
		else
			histData.prctileConeIsolatingBPIindex(iBin,:) = nan(1,numel(thePercentileValues));
		end

		histData.meanConeIsolatingBPIindex(iBin) = mean(d(:));
	end

	bpiMap = bpiMap / max(bpiMap(:));
	bpiMap(bpiMap>maxVisualizedBPIDensity) = 1.0;

	if (centerConeDominance == cMosaic.LCONE_ID)
		centerColor = RGCMosaicConstructor.constants.LcenterColor;
		centerColorLeeShapley = 'none';
		if (superimposeLeeShapleyData)
			dataOut = RGCMosaicConstructor.publicData.LeeShapley.midgetBandPassIndices();
			LeeShapleyBPIdata = dataOut('L-center');
		end
	elseif (centerConeDominance == cMosaic.MCONE_ID)
		centerColor = RGCMosaicConstructor.constants.McenterColor;
		centerColorLeeShapley = 'none';
		if (superimposeLeeShapleyData)
			dataOut = RGCMosaicConstructor.publicData.LeeShapley.midgetBandPassIndices();
			LeeShapleyBPIdata = dataOut('M-center');
		end
	else
		centerColor = RGCMosaicConstructor.constants.achromaticColor;
		centerColorLeeShapley = 'none';
		if (superimposeLeeShapleyData)
			dataOut = RGCMosaicConstructor.publicData.LeeShapley.midgetBandPassIndices();
			LeeShapleyBPIdataLcenterData = dataOut('L-center');
			LeeShapleyBPIdataMcenterData = dataOut('M-center');
			LeeShapleyBPIdata.x = cat(1, LeeShapleyBPIdataLcenterData.x(:),  LeeShapleyBPIdataMcenterData.x(:));
			LeeShapleyBPIdata.y = cat(1, LeeShapleyBPIdataLcenterData.y(:),  LeeShapleyBPIdataMcenterData.y(:));
		end
	end
	achromaticColor = RGCMosaicConstructor.constants.achromaticColor;

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = 'validation';


	renderDensityPlot = ~isempty(pdfFileName4);
	renderPercentilesPlot = ~isempty(pdfFileName5);
	renderScatterPlot = ~isempty(pdfFileName);
	renderHistogramPlot = ~isempty(pdfFileName2);
	renderSurroundConePurityVsBPIplot = ~isempty(pdfFileName3);

	if (renderDensityPlot)
		% The density plot of acrom vs cone isolating BPIs
		ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
		hFig = figure(figNo); clf;
	    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	    ax = theAxes{1,1};

	    bpiLevels = 0.05: 0.1: 1.0;
	    %imagesc(ax, bpiBins+bpiInterval, bpiBins+bpiInterval, bpiMap);

	    [X,Y] = meshgrid(bpiBins+bpiInterval, bpiBins+bpiInterval);
	    [Xq,Yq] = meshgrid(0.0:0.01:1, 0.0:0.01:1);
	    Vq = interp2(X,Y,bpiMap,Xq,Yq, 'pchip');

	    contourf(ax, Xq,Yq, Vq, bpiLevels, 'LineWidth', 1.0, 'Color', [0 0 0.0]);
	    hold(ax, 'on');

	    if (superimposeLeeShapleyData)
	    	scatter(ax, LeeShapleyBPIdata.x, LeeShapleyBPIdata.y, superimposeLeeShapleyDataMarkerSize^2, ...
	    		'Marker', 's', ...
	    		'LineWidth', 0.75, ...
	    		'MarkerFaceColor', centerColor, ...
	    		'MarkerEdgeColor', [0 0 0], ...
	    		'MarkerFaceAlpha', 0.6);
	    end

	    axis (ax, 'square');
	    axis(ax, 'xy');
	    XLims = [0 1.0]; YLims = [0.3 1.0];
		grid(ax, 'on');
		set(ax, 'CLim', [0 1.1]);
		set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
			    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
	    xlabel(ax,'achromatic gratings');
	    ylabel(ax,'cone-isolating gratings');
	    cMap = brewermap(512, 'greys');
	    cMap(1,:) = [1 1 1];

	    colormap(ax, cMap);

	    % Finalize figure using the Publication-Ready format
	    PublicationReadyPlotLib.applyFormat(ax,ff);
	    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName4);
	    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end


	if (renderPercentilesPlot)
	    % The prctile plot of acrom vs cone isolating BPIs
		ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
		hFig = figure(figNo+1); clf;
	    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	    ax = theAxes{1,1};

	    hold(ax, 'on');

		if (superimposeLeeShapleyData)	
	    	p6 = scatter(ax, LeeShapleyBPIdata.x, LeeShapleyBPIdata.y, superimposeLeeShapleyDataMarkerSize^2, ...
	    		'Marker', 'o', ...
	    		'LineWidth', 1.0, ...
	    		'MarkerFaceColor', centerColor.^2, ...
	    		'MarkerEdgeColor', [1 1 1], ...
	    		'MarkerFaceAlpha', 0.7);
	    end

	    if (~onlyDepictLeeShapleyData)

	    	newStyleDistributionVisualization = true;
	    	if (newStyleDistributionVisualization)
	    		for iP = 1:size(histData.prctileConeIsolatingBPIindex,2)/2-1
	    			xx = histData.achromaticBPIs + histData.achromaticBPIs * bpiInterval;
	    			zLevel1data.x = xx(:);
	    			zLevel2data.x = flipud(xx(:));
	    			zLevel1data.y = histData.prctileConeIsolatingBPIindex(:,iP);
	    			zLevel2data.y = histData.prctileConeIsolatingBPIindex(:,size(histData.prctileConeIsolatingBPIindex,2)-iP);
		    		X  = cat(1, zLevel1data.x, zLevel2data.x);            
					Yp = cat(1, zLevel1data.y, flipud(zLevel2data.y));
					validIndices = find(~isnan(Yp));
					Yp = Yp(validIndices);
					X = X(validIndices);
					v = [X(:) Yp(:)];
					f = 1:numel(X);
					patch(ax, 'Faces',f,'Vertices',v, 'FaceColor', [50 100 170]/255, 'FaceAlpha',0.15, 'EdgeColor', [0 0 1], 'EdgeAlpha', 0.2);
				end

	    	else
				% The 5%
				plot(ax, histData.achromaticBPIs, histData.prctileConeIsolatingBPIindex(:,1), ':', 'LineWidth', 2, ...
					'MarkerSize', markerSize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0 0], 'Color', [1 0 0]);

				% The 25%
				plot(ax, histData.achromaticBPIs, histData.prctileConeIsolatingBPIindex(:,2), '--', 'LineWidth', 2, ...
					'MarkerSize', markerSize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0 0], 'Color', [1 0 0]);

				% The median
				plot(ax, histData.achromaticBPIs, histData.prctileConeIsolatingBPIindex(:,3), '-', 'LineWidth', 3, ...
					'MarkerSize', markerSize, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'Color', [1 0 0]);

				% The 75%
				plot(ax, histData.achromaticBPIs, histData.prctileConeIsolatingBPIindex(:,4), '--', 'LineWidth', 2, ...
					'MarkerSize', markerSize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0 0], 'Color', [1 0 0]);

				% The 95%
				plot(ax, histData.achromaticBPIs, histData.prctileConeIsolatingBPIindex(:,5), ':', 'LineWidth', 2, ...
					'MarkerSize', markerSize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0 0], 'Color', [1 0 0]);
			end

		end

		if (superimposeLeeShapleyData)	
	    	p6 = scatter(ax, LeeShapleyBPIdata.x, LeeShapleyBPIdata.y, superimposeLeeShapleyDataMarkerSize^2, ...
	    		'Marker', 'o', ...
	    		'LineWidth', 1.0, ...
	    		'MarkerFaceColor', centerColor.^2, ...
	    		'MarkerEdgeColor', [1 1 1], ...
	    		'MarkerFaceAlpha', 0.4);
	    end

		if (showLegends)
		    if (onlyDepictLeeShapleyData)
		    	theLegends{1} = 'macaque';
			    legend(ax, p6, theLegends, 'Location', 'SouthEast');
		    else
			    if (superimposeLeeShapleyData)
			    	theLegends = {};
			    	
			    	for i = 1:numel(thePercentileValues)
			    		theLegends{numel(theLegends)+1} = sprintf('synthetic (%2.0f%%)', thePercentileValues(i));
					end
					theLegends{numel(theLegends)+1} = 'macaque';
			    	legend(ax, [p1 p2 p3 p4 p5 p6], theLegends, 'Location', 'SouthEast');
			    end
			end
		end


	    axis (ax, 'square');
	    XLims = [0 1.0]; YLims = [0 1.0];
		grid(ax, 'on');
		set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
			    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
	    xlabel(ax,'achromatic gratings');
	    ylabel(ax,'center one-isolating gratings');


	    % Finalize figure using the Publication-Ready format
	    PublicationReadyPlotLib.applyFormat(ax,ff);
	    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName5);
	    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end



	if (renderScatterPlot)
		% The scatter plot of acrom vs cone isolating BPIs
		ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
		hFig = figure(figNo+1); clf;
	    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	    ax = theAxes{1,1};

	    scatter(ax, achromaticBPIs, coneIsolatingBPIs, markerSize^2, ...
	    	'MarkerFaceColor', centerColor, 'MarkerEdgeColor', centerColor, ...
	    	'MarkerFaceAlpha', markerFaceAlpha, ...
	    	'MarkerEdgeAlpha', markerEdgeAlpha, ...
	    	'LineWidth', markerLineWidth);

	    hold(ax, 'on');

	    if (superimposeLeeShapleyData)	
	    	scatter(ax, LeeShapleyBPIdata.x, LeeShapleyBPIdata.y, superimposeLeeShapleyDataMarkerSize^2, ...
	    		'Marker', 's', ...
	    		'MarkerFaceColor', [0.7 0.7 0.7], ...
	    		'MarkerEdgeColor', [0 0 0], ...
	    		'MarkerFaceAlpha', 0.5);
	    end

	    if (superimposeLeeShapleyData)
	    	legend(ax, {'synthetic', 'macaque'}, 'Location', 'SouthWest');
	    end

	    axis (ax, 'square');
	    XLims = [0 1.0]; YLims = [0 1.0];
		grid(ax, 'on');
		set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
			    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
		set(ax, 'CLim', [0 1.2]);
	    xlabel(ax,'achromatic gratings');
	    ylabel(ax,'cone-isolating gratings');

	    % Finalize figure using the Publication-Ready format
	    PublicationReadyPlotLib.applyFormat(ax,ff);
	    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, strrep(pdfFileName, '.pdf', '.png'));
	    NicePlot.exportFigToPNG(thePDFfileName, hFig,  300);
	end



	if (renderHistogramPlot)
	    % The marginal histograms
	    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
		hFig = figure(figNo+1); clf;
	    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	    ax = theAxes{1,1};

		bpiBins = 0:0.02:1.0;
	    faceAlpha = 0.6;
	    lineWidth = 1.5;

	    bpiBins = 0:0.1:1.0;

	    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, coneIsolatingBPIs, bpiBins, ...
	    	 centerColor, centerColor*0.5, faceAlpha, lineWidth, ...
	    	 'relativeHeight', 1.0, ...
	    	 'binWidth', 0.7);
	    hold(ax, 'on');

	    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, achromaticBPIs, bpiBins, ...
	    	 achromaticColor, achromaticColor*0.5, faceAlpha, lineWidth, ...
	    	 'relativeHeight', 1.0, ...
	    	 'binWidth', 0.9);

	    if (superimposeLeeShapleyData)
			dataOut = RGCMosaicConstructor.publicData.LeeShapley.midgetBandPassIndices();
			LeeShapleyBPIdataLcenters = dataOut('L-center');
			LeeShapleyBPIdataMcenters = dataOut('M-center');
			LeeShapleyAchromaticBPIs = cat(1,LeeShapleyBPIdataLcenters.x(:), LeeShapleyBPIdataMcenters.x(:));
			LeeShapleyConeIsolatingBPIs = cat(1,LeeShapleyBPIdataLcenters.y(:), LeeShapleyBPIdataMcenters.y(:));

		    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, LeeShapleyAchromaticBPIs, bpiBins, ...
		    	 achromaticColor*0.5, [0 0 0], faceAlpha, lineWidth, ...
		    	 'relativeHeight', 1.0, ...
		    	 'binWidth', 0.7, ...
		    	 'plotWithNegativePolarity', true);

		    RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, LeeShapleyConeIsolatingBPIs, bpiBins, ...
		    	 centerColor*0.5, [0 0 0], faceAlpha, lineWidth, ...
		    	 'relativeHeight', 1.0, ...
		    	 'binWidth', 0.9, ...
		    	 'plotWithNegativePolarity', true);

	    	legend(ax, {...
	    		'cone isolating (synthetic)', ...
	    		'achromatic (synthetic)', ...
	    		'cone isolating (macaque)', ...
	    		'achromatic (macaque)'}, ...
	    		'Location', 'SouthWest', 'FontSize', 12);
	    end

		XLims = [0 1];
		YLims = get(ax, 'YLim');
		if (superimposeLeeShapleyData)
			YLims = max(abs(YLims(:)))*[-1 1];
		end

		axis(ax, 'square');
		set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
			'YLim', YLims, 'YTickLabel', {});
	    xlabel(ax,'BPI');
	    ylabel(ax, 'frequency')
	    % Finalize figure using the Publication-Ready format
	    PublicationReadyPlotLib.applyFormat(ax,ff);
	    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName2);
	    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end

	if (renderSurroundConePurityVsBPIplot)
	    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
		hFig = figure(figNo+2); clf;
	    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	    ax = theAxes{1,1};

		scatter(ax, surroundConePurities, coneIsolatingBPIs, 100, ...
	    	'MarkerFaceColor', centerColor, 'MarkerEdgeColor', centerColor*0.5, 'MarkerFaceAlpha', 0.3);

		axis (ax, 'square');
	    XLims = [0 1.0]; YLims = [0 1.0];
		grid(ax, 'on');
		set(ax, 'XLim', XLims, 'XTick', [0 1/6 1/3 1/2 2/3 5/6 1], 'XTickLabel', {'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'}, ...
			    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
	    xlabel(ax,'surround cone purity');
	    ylabel(ax,'BPI (cone-isolating gratings)');

	    % Finalize figure using the Publication-Ready format
	    PublicationReadyPlotLib.applyFormat(ax,ff);
	    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName3);
	    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end


end
