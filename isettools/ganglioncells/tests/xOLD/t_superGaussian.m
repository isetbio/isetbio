function t_superGaussian()

	stDev = 0.2;

	exponentsVisualized = [2 10];
	cLUT = [...
		1 0 0; ...
	    0.2 0.2 0.2 ...
	    ];

	singleGaussian(1, exponentsVisualized, stDev, cLUT, 'poolingWeights.pdf');

	for i = 1:numel(exponentsVisualized)
		neighboringGaussians(1+i, exponentsVisualized(i), stDev, cLUT(i,:), sprintf('RFcenterOverlaps_Exp_%2.2f.pdf', exponentsVisualized(i)));
	end
end

function neighboringGaussians(figNo, exponentVisualized, stDev, faceColor, pdfFileName)
	% Generate figure
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    retinalSpace = -1:0.001:1;

    superGaussianLeft = exp(-0.5*(abs((retinalSpace-2*stDev)/stDev)).^exponentVisualized);
    superGaussianCenter = exp(-0.5*(abs((retinalSpace)/stDev)).^exponentVisualized);
    superGaussianRight = exp(-0.5*(abs((retinalSpace+2*stDev)/stDev)).^exponentVisualized);

    superGaussianPosSigma = exp(-0.5*(abs(stDev/stDev)).^10);
    superGaussianNegSigma = exp(-0.5*(abs(-stDev/stDev)).^10);

    edgeColor = faceColor;
	faceColor = 0.5*faceColor + [0.5 0.5 0.5]
    
    faceAlpha = 0.5;
    lineWidth = 2.0;
    lineStyle = ':';
    
    lineWidth = 1.0;
    lineStyle = '-';
    RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(theAxes{1,1}, ...
    	retinalSpace, superGaussianLeft, superGaussianLeft*0, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);

    hold(theAxes{1,1}, 'on')

    lineWidth = 1.0;
    lineStyle = '-';
    RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(theAxes{1,1}, ...
    	retinalSpace, superGaussianRight, superGaussianRight*0, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);

    lineWidth = 2.0;
    lineStyle = '-';
    RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(theAxes{1,1}, ...
    	retinalSpace, superGaussianCenter, superGaussianCenter*0, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);

    plot(theAxes{1,1}, stDev, superGaussianPosSigma , 'ro', 'LineWidth', 1.5, ...
    	'MarkerSize', 16, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    plot(theAxes{1,1}, -stDev, superGaussianNegSigma , 'ro', 'LineWidth', 1.5, ...
    	'MarkerSize', 16, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);

    % Axes scaling and ticks
    set(theAxes{1,1}, 'XTick', -1:0.5:1);
    set(theAxes{1,1}, 'YTick', 0:0.2:1.0);
    
    % Finalize figure using the Publication-Ready format
    XLims = [-1 1]; YLims = [0 1.05];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'retinal space (degs)', 'cone pooling weight');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
end

function singleGaussian(figNo, exponentIvisualized, stDev, cLUT, pdfFileName)
	% Generate figure
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    retinalSpace = -1:0.001:1;

    for i = 1:numel(exponentIvisualized)
    	superGaussian(i,:) = exp(-0.5*(abs(retinalSpace/stDev)).^exponentIvisualized(i) );
    end

    superGaussianPosSigma = exp(-0.5*(abs(stDev/stDev)).^10);
    superGaussianNegSigma = exp(-0.5*(abs(-stDev/stDev)).^10);

    for exponentI = 1:size(superGaussian,1)
    	plot(theAxes{1,1}, retinalSpace, superGaussian(exponentI,:), ...
    		'LineWidth', 3.0, 'Color', cLUT(exponentI,:));
    	hold(theAxes{1,1}, 'on')
    end

    plot(theAxes{1,1}, stDev, superGaussianPosSigma , 'ro', 'LineWidth', 1.5, ...
    	'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
    plot(theAxes{1,1}, -stDev, superGaussianNegSigma , 'ro', 'LineWidth', 1.5, ...
    	'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);

    % Axes scaling and ticks
    set(theAxes{1,1}, 'XTick', -1:0.5:1);
    set(theAxes{1,1}, 'YTick', 0:0.2:1.0);
    
    % Finalize figure using the Publication-Ready format
    XLims = [-1 1]; YLims = [0 1.05];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'retinal space (degs)', 'cone pooling weight');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
end
