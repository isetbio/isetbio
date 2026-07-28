function compareHumanToMacaque

    degsRange = 120*[-1 1];
    degsTicks = -120:30:120;
    MMsRange = 30*[-1 1];
    MMsTicks = -30:5:30;

    xRangeMMs = MMsRange(1):0.1:MMsRange(2);
    xRangeDegsHuman = RGCmodels.Watson.convert.rhoMMsToDegs(xRangeMMs);
    xRangeDegsMacaque = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', xRangeMMs);


    micronsPerDegHumanRetina = 1e3*xRangeMMs./xRangeDegsHuman;
    micronsPerDegMacaqueRetina = 1e3*xRangeMMs./xRangeDegsMacaque;
    maxMicronsPerDegMacaqueRetina = max(micronsPerDegMacaqueRetina);
    maxMicronsPerDegHumanRetina = max(micronsPerDegHumanRetina);
    

    figNo = 1;
    ff = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
    ff.axisFontSize = 24;
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
    p1 = plot(ax, xRangeMMs, xRangeDegsMacaque, 'r-', 'LineWidth', 1.5); 
    hold(ax, 'on');
    p2 = plot(ax, xRangeMMs, xRangeDegsHuman, 'b-', 'LineWidth', 1.5);
    %plot(ax, xRangeMMs, xRangeMMs * 1e3 / maxMicronsPerDegHumanRetina, 'b--', 'LineWidth', 1.5);
    set(ax, 'XLim', [0 max(MMsRange)], 'XTick', MMsTicks, 'YLim', [0 max(degsRange)], 'YTick', degsTicks);
    xlabel(ax, 'retinal eccentricity (mm)');
    ylabel(ax,'retinal eccentricity (degs)');
    legend(ax, [p1 p2], {'macaque', 'human'}, 'Location', 'SouthEast');

    % Finalize figure using the Publication-Ready format
	ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, xLims, yLims);
    NicePlot.exportFigToPDF('DegsHumanVsMacaque.pdf', hFig, 300);



    figNo = 2;
    ff = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
    ff.axisFontSize = 24;
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
    p1 = plot(ax, xRangeMMs, micronsPerDegMacaqueRetina, 'r-', 'LineWidth', 1.5); 
    hold(ax, 'on');
    p2 = plot(ax, xRangeMMs, micronsPerDegHumanRetina, 'b-', 'LineWidth', 1.5);
    set(ax, 'XLim', [0 max(MMsRange)], 'XTick', MMsTicks, 'YLim', [0 300], 'YTick', 0:50:300);
    xlabel(ax, 'retinal eccentricity (mm)');
    ylabel(ax,'microns/deg');
    legend(ax, [p1 p2], {'macaque', 'human'}, 'Location', 'SouthEast');

    % Finalize figure using the Publication-Ready format
	ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, xLims, yLims);
    NicePlot.exportFigToPDF('DegsHumanVsMacaque2.pdf', hFig, 300);


    figNo = 3;
    ff = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
    ff.axisFontSize = 24;
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
    plot(ax, xRangeDegsHuman, xRangeDegsMacaque, 'k-', 'LineWidth', 3.0); 
    hold(ax, 'on')
    plot(ax, xRangeDegsHuman, xRangeDegsMacaque, 'c-', 'LineWidth', 1.5);
    plot(ax, xRangeDegsHuman, xRangeDegsHuman, 'k--', 'LineWidth', 1.0);
    axis(ax, 'square');
    set(ax, 'XLim', [0 max(degsRange)], 'XTick', degsTicks, 'YLim', [0 max(degsRange)], 'YTick', degsTicks);
    xlabel(ax, 'human retinal eccentricity (degs)');
    ylabel(ax, 'macaque retinal eccentricity (degs)');

    % Finalize figure using the Publication-Ready format
	ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, xLims, yLims);
    NicePlot.exportFigToPDF('DegsHumanVsMacaque3.pdf', hFig, 300);

    
end
