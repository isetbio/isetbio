%
% RGCMosaicAnalyzer.visualize.mosaicRetinalDiametersAgainstMacaqueInVitroRFdiameters(..)
%

function mosaicRetinalDiametersAgainstMacaqueInVitroRFdiameters(figNo, ff,...
	radialEccMMs, minorSigmaMicrons, majorSigmaMicrons, sigmaMicrons, ...
    displayMinorAndMajorAxisDiameters, employLogYaxis, employLogXaxis, ...
    pdfFileName)

	% Load the EJ data
    [dataOut, sourceDocument] = RGCMosaicConstructor.publicData.ChichilniskyLab.rfDiameters();
    d85mm = dataOut('8.5mm');
    d35mm = dataOut('3.5mm');

    % Plot the RF diameters together with the EJ data
    hFig = figure(figNo); clf;
    ff.backgroundColor = [1 1 1];
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
    hold(ax, 'on');

    if (displayMinorAndMajorAxisDiameters)
    	pMinor = scatter(ax, radialEccMMs, 2*minorSigmaMicrons, 81, 'MarkerFaceColor', [0.95 0.3 0.8], 'MarkerFaceAlpha', 0.3, ...
        	'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.5,'MarkerEdgeAlpha', 0.2);
    	pMajor = scatter(ax, radialEccMMs, 2*majorSigmaMicrons, 81, 'MarkerFaceColor', [0 0.5 1], 'MarkerFaceAlpha', 0.3, ...
        	'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.2);
    end

    darkGreenColor = [0.0 0.3 0];
    greenColor = [0.0 1.0 0.2];
    blueColor = [0.0 0.4 1.0];
    lightGrayColor = [0.85 0.85 0.85];
    darkGrayColor = [0.0 0.0 0.0];
    pSynthetic = scatter(ax, radialEccMMs, 2*sigmaMicrons, 4, 'MarkerFaceColor', darkGrayColor, 'MarkerFaceAlpha', 1, ...
        'MarkerEdgeColor', 'none', 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.9);

    scatter(ax, 3.5+0.2*randn(size(d35mm.rfDiameterMicrons)), d35mm.rfDiameterMicrons, 121, 'Marker', 's', ...
        'MarkerFaceColor', greenColor, 'MarkerFaceAlpha', 0.4, 'LineWidth', 0.75, ...
        'MarkerEdgeColor', darkGreenColor, 'MarkerEdgeAlpha', 1);

    pInVitro = scatter(ax, 8.5+0.2*randn(size(d85mm.rfDiameterMicrons)), d85mm.rfDiameterMicrons, 121, 'Marker', 's', ...
        'MarkerFaceColor', greenColor, 'MarkerFaceAlpha', 0.4, 'LineWidth', 0.75, ...
        'MarkerEdgeColor', darkGreenColor,  'MarkerEdgeAlpha', 1);


    if (displayMinorAndMajorAxisDiameters)
    	legend(ax, [pMinor pMajor pSynthetic pInVitro], {'synthetic mRGCs (minor)', 'synthetic mRGCs (major)', 'synthetic mRGCs', sprintf('macaque mRGCs, %s',sourceDocument(1:24))}, ...
    		'Location', 'NorthWest', 'FontSize', 12, 'NumColumns', 1);
    else
    	legend(ax, [pSynthetic pInVitro], {'synthetic mRGCs', sprintf('macaque mRGCs, %s',sourceDocument(1:24))}, ...
    		'Location', 'NorthWest', 'FontSize', 12, 'NumColumns', 1);
    end

    if (employLogXaxis)
    	set(ax, 'XLim', [0.01 10], 'XTick', [0.01 0.03 0.1 0.3 1 3 10], 'XScale', 'log');
    else
    	set(ax, 'XLim', [0 10], 'XTick', 0:1:10, 'XScale', 'linear');
    end

    if (employLogYaxis)
    	set(ax, 'YScale', 'log', 'YLim', [1 80], 'YTick', [1 2 4 8 16 32 64]);
    else
    	set(ax, 'YScale', 'linear', 'YLim', [0 80], 'YTick', 0:10:100);
    end
    
    xlabel(ax, 'radial eccentricity (mm)');
    ylabel(ax, 'RF diameter (microns)');

    ff.box = 'off';
    ff.tickDir = 'in';
    ff.grid = 'on';
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    % axis(ax, 'square');

    % Export the PDF
    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);

    % Using the export_fig which supports transparent background
    %export_fig(hFig, thePDFfileName, '-pdf', '-native');
end
