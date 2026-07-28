%
% RGCMosaicAnalyzer.visualize.correspondenceScatterPlot(..)
%


function correspondenceScatterPlot(figNo, ...
        theVisualSpaceReferredData, theRetinalSpaceReferredData, ...
        markerSize, markerType, markerColor, markerFaceAlpha, markerEdgeAlpha, ...
        theRange, theTicks, theMeasurementString, pdfFileName)

    % The density plot of acrom vs cone isolating BPIs
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    scatter(ax, theVisualSpaceReferredData, theRetinalSpaceReferredData, markerSize^2, ...
                'Marker', markerType, ...
                'LineWidth', ff.lineWidth/2, ...
                'MarkerFaceColor', markerColor, ...
                'MarkerEdgeColor', markerColor*0.5, ...
                'MarkerFaceAlpha', markerFaceAlpha, ...
                'MarkerEdgeAlpha', markerEdgeAlpha);

    axis(ax, 'square');
    set(ax, 'XLim', theRange, 'XTick', theTicks, 'YLim', theRange, 'YTick', theTicks);
    xlabel(ax, sprintf('%s (visual space - referred)',theMeasurementString));
    ylabel(ax, sprintf('%s (retinal space - referred)',theMeasurementString));

    % Finalize figure using the Publication-Ready format
    ff.legendBox = 'on';
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, theRange, theRange);

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);

end

