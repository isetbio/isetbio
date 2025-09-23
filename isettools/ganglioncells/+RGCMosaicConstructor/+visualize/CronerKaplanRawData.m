function CronerKaplanRawData(varargin)

    p = inputParser;
    p.addParameter('surroundToCenterRcRatioTolerance', 15, @(x)(isempty(x)||isscalar(x)));
    p.addParameter('surroundToCenterIntegratedSensitivityRatioTolerance', 30, @(x)(isempty(x)||isscalar(x)));
    p.addParameter('surroundToCenterRcRatioMultiplier', 1.0, @(x)(isempty(x)||isscalar(x)));
    p.addParameter('surroundToCenterIntegratedSensitivityRatioMultiplier', 1.0, @(x)(isempty(x)||isscalar(x)));
    p.parse(varargin{:});
    surroundToCenterRcRatioTolerance = p.Results.surroundToCenterRcRatioTolerance;
    surroundToCenterIntegratedSensitivityRatioTolerance = p.Results.surroundToCenterIntegratedSensitivityRatioTolerance;
    surroundToCenterRcRatioMultiplier = p.Results.surroundToCenterRcRatioMultiplier;
    surroundToCenterIntegratedSensitivityRatioMultiplier = p.Results.surroundToCenterIntegratedSensitivityRatioMultiplier;

	% RGCMosaicConstructor.visualize.CronerKaplanRawData
	typicalSTF = RGCmodels.CronerKaplan.digitizedData.typicalRGCSTF();
	%plot(typicalSTF.sfCPD, typicalSTF.amplitude, 'k-');

	[eccentricityDegs, RcRsRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    RsRcRatios = 1./RcRsRatios;
    theTrend = prctile(RsRcRatios, 50) +  0*eccentricityDegs;
    theDetrendedRatios = RsRcRatios - theTrend;

    theTrend = theTrend * surroundToCenterRcRatioMultiplier;
    lowRangeRatios = theTrend + prctile(theDetrendedRatios, 50-surroundToCenterRcRatioTolerance/2);
    highRangeRatios = theTrend + prctile(theDetrendedRatios, 50+surroundToCenterRcRatioTolerance/2);

    [~,idx] = sort(eccentricityDegs, 'ascend');
    eccentricityDegs = eccentricityDegs(idx);
    RsRcRatios = RsRcRatios(idx);
    theTrend = theTrend(idx);
    lowRangeRatios = lowRangeRatios(idx);
    highRangeRatios  = highRangeRatios (idx);

    figNo = 1; thePDFfileName = 'RsRcRatio.pdf';
    [hFig, ax, ff] = startThePlot(figNo);
    XLims = [0 40];
    YLims = [1 20];
    
    xx = eccentricityDegs(:);
    yy = lowRangeRatios(:);
    xx = cat(1, eccentricityDegs(:), eccentricityDegs(end));
    yy = cat(1, yy(:), highRangeRatios(end));
    xx = cat(1, xx(:), flipud(eccentricityDegs(:)));
    yy = cat(1, yy(:), flipud(highRangeRatios(:)));
    zz = -10*eps*ones(size(yy));
    c = 0.5*([0.0 114 189]/255+ [0.5 0.5 0.5]);
    patch(ax,xx,yy,zz, 'FaceColor', c, 'EdgeColor', c*0.5, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
    hold(ax, 'on');

    
    plot(ax, eccentricityDegs, RsRcRatios, 'ks', 'MarkerSize', ff.markerSize, ...
    	'MarkerEdgeColor', [0.0 114 189]/255*0.5, 'MarkerFaceColor', [0.0 114 189]/255, 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    plot(ax, eccentricityDegs, theTrend, 'k--', 'LineWidth', ff.lineWidth*2, 'Color', [0.0 114 189]/255*0.7);
    set(ax, 'XTick', 0:10:100, 'YTick', 0:5:20);
    hold(ax, 'off');
    set(ax, 'XLim', XLims, 'YLim', YLims);
    xlabel(ax, 'temporal equivalent eccentricity (degs)');
    ylabel(ax, 'Rs/Rc ratio');
    axis(ax, 'square');
    finishThePlot(ax, ff, hFig, XLims, YLims, thePDFfileName);


    [eccentricityDegs, SCintegratedSensitivityRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
    theTrend = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(eccentricityDegs);
    theDetrendedRatios = SCintegratedSensitivityRatios - theTrend;

    theTrend = theTrend * surroundToCenterIntegratedSensitivityRatioMultiplier;
    lowRangeRatios = theTrend + prctile(theDetrendedRatios, 50-surroundToCenterIntegratedSensitivityRatioTolerance/2);
    highRangeRatios = theTrend + prctile(theDetrendedRatios, 50+surroundToCenterIntegratedSensitivityRatioTolerance/2);

    [~,idx] = sort(eccentricityDegs, 'ascend');
    eccentricityDegs = eccentricityDegs(idx);
    SCintegratedSensitivityRatios = SCintegratedSensitivityRatios(idx);
    theTrend = theTrend(idx);
    lowRangeRatios = lowRangeRatios(idx);
    highRangeRatios  = highRangeRatios (idx);

    figNo = 2; thePDFfileName = 'intSCRatio.pdf';
    [hFig, ax, ff] = startThePlot(figNo);
    XLims = [0 40];
    YLims = [0 1];
    

    xx = eccentricityDegs(:);
    yy = lowRangeRatios(:);
    xx = cat(1, eccentricityDegs(:), eccentricityDegs(end));
    yy = cat(1, yy(:), highRangeRatios(end));
    xx = cat(1, xx(:), flipud(eccentricityDegs(:)));
    yy = cat(1, yy(:), flipud(highRangeRatios(:)));
    zz = -10*eps*ones(size(yy));
    c = 0.5*([218 83 15]/255+ [0.5 0.5 0.5]);
    patch(ax,xx,yy,zz, 'FaceColor', c, 'EdgeColor', c*0.5, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
    hold(ax, 'on');

    plot(ax, eccentricityDegs, SCintegratedSensitivityRatios, 'ks', 'MarkerSize', ff.markerSize, 'MarkerEdgeColor', [218 83 15]/255*0.5, 'MarkerFaceColor', [218 83 15]/255, 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    
    plot(ax, eccentricityDegs, theTrend , '--', 'Color', [218 83 15]/255*0.7, 'LineWidth', ff.lineWidth*2, 'Color', [1 0.3 0.0]*0.5);
    set(ax, 'XTick', 0:10:100, 'YTick', 0:0.2:1);
    hold(ax, 'off');
    set(ax, 'XLim', XLims, 'YLim', YLims);
    xlabel(ax, 'temporal equivalent eccentricity (degs)');
    ylabel(ax, 'intS/intC ratio');
    axis(ax, 'square');
    finishThePlot(ax, ff, hFig, XLims, YLims, thePDFfileName);
end

function [hFig, ax, ff] = startThePlot(figNo)
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
end

function finishThePlot(ax,ff, hFig, XLims, YLims, thePDFfileName)
	% Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = '';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFfileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end
