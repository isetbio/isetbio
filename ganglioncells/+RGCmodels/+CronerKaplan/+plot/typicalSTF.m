% Generate a typical STF of a macaque RGC

function [hFig, hFig2, hFig3] = typicalSTF(figNo)
	typicalSTFdata = RGCmodels.CronerKaplan.digitizedData.typicalRGCSTF();

	% Generate figure
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    % Fit DoG model to the data
    radialTemporalEquivalentEccentricityDegs = 5;
    RsDegs = RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialTemporalEquivalentEccentricityDegs);
	RcDegsInitial = RsDegs/RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

	RcDegs = struct(...
		'initial', RcDegsInitial, ...
    	'low', RcDegsInitial/100, ...
    	'high', RcDegsInitial*100);

    initialKc = RGCmodels.CronerKaplan.constants.centerPeakSensitivityFromCharacteristicRadiusDegsForPcells(RcDegs.initial);
	Kc = struct(...
				'initial', initialKc, ... 
	        	'low', initialKc/1000, ...
	        	'high', initialKc*1000);

	intStoCsens = struct(...
		'initial', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccentricityDegs), ...
		'low', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(0)*0.2, ...
		'high', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(30)*5);


	RsToRc = struct(...
		'initial', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio, ...
		'low', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio/10, ...
		'high', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio*10);

    DoGparams.initialValues = [Kc.initial   intStoCsens.initial    RsToRc.initial    RcDegs.initial];
	DoGparams.lowerBounds   = [Kc.low       intStoCsens.low        RsToRc.low        RcDegs.low];
	DoGparams.upperBounds   = [Kc.high      intStoCsens.high       RsToRc.high       RcDegs.high];
	DoGparams.names         = {'Kc',        'intStoCsens',         'RsToRc',         'RcDegs'};
	DoGparams.scaling       = {'log',       'linear',           'linear',         'linear'};

	axDoGFitToCompositeSTF = [];
	normFactor = 1;
	typicalSTFdata.amplitude = typicalSTFdata.amplitude * 0.9;
	multiStartsNum = 32;
    [DoGparams, theFittedSTF] = RGCMosaicConstructor.helper.fit.genericDifferentOfGaussiansToCompositeSTF(...
	    		DoGparams, typicalSTFdata.sfCPD, typicalSTFdata.amplitude, ...
	    		axDoGFitToCompositeSTF, normFactor, multiStartsNum);

    
    idx = find(strcmp(DoGparams.names,'Kc'));
	KcBest = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'intStoCsens'));
	intStoCsensBest = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'RsToRc'));
	RsToRcBest = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'RcDegs'));
	RcDegsBest = DoGparams.finalValues(idx);

	plot(theAxes{1,1}, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, '-', ...
		'Color', [0 0 0], 'LineWidth', ff.lineWidth*2);
	hold(theAxes{1,1}, 'on');
	plot(theAxes{1,1}, theFittedSTF.sfHiRes, theFittedSTF.centerSTFHiRes, 'k-', 'LineWidth', ff.lineWidth);
	plot(theAxes{1,1}, theFittedSTF.sfHiRes, theFittedSTF.surroundSTFHiRes, 'k--', 'LineWidth', ff.lineWidth);
	plot(theAxes{1,1}, typicalSTFdata.sfCPD, typicalSTFdata.amplitude, 'o', ...
		'MarkerSize', ff.markerSize*1.3, ...
		'MarkerFaceColor', [0.5 0.85 1], 'MarkerEdgeColor', 0.5*[0.5 0.85 1], 'LineWidth', ff.lineWidth);

	% Axes scaling and ticks
    set(theAxes{1,1}, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(theAxes{1,1}, 'YScale', 'linear', ...
        'YTick', 0:0.2:1, ...
        'YTickLabel', 0:0.2:1);

    % Finalize figure using the Publication-Ready format
    XLims = [0.01 30]; YLims = [0 1];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'spatial frequency (c/deg)', 'response (normalized)');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

    figNo = figNo+1;
    hFigTmp = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFigTmp,ff);
    ax = theAxes{1,1};
    RGCMosaicConstructor.visualize.fittedModelParams(ax, DoGparams, 'DoG');


    [eccentricityDegsCKdata, RcToRsCKdata] = ...
		RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
	RsRcCKdata = 1./RcToRsCKdata;

    figNo = figNo+1;
    hFig2 = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig2,ff);
    ax = theAxes{1,1};

    plot(theAxes{1,1}, eccentricityDegsCKdata, RsRcCKdata, 'ks', ...
		'MarkerSize', ff.markerSize, ...
		'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'LineWidth', ff.lineWidth*0.75);

    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1}, 10, RsToRcBest, 'o', 'MarkerEdgeColor', 0.5*[0.5 0.85 1], ...
    	'MarkerFaceColor', [0.5 0.85 1], 'MarkerSize', ff.markerSize*1.25, 'LineWidth', 2);
    set(theAxes{1,1}, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);

	% Finalize figure using the Publication-Ready format
    XLims = [0.1 40]; YLims = [0 20];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'eccentricity (degs)', 'Rs/Rc');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);


    [eccentricityDegsCKdata, intSCsensCKdata] = ...
		RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
 
 	figNo = figNo+1;
    hFig3 = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig3,ff);
    ax = theAxes{1,1};

    plot(theAxes{1,1}, eccentricityDegsCKdata, intSCsensCKdata, 'ks', ...
		'MarkerSize', ff.markerSize, ...
		'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0 0 0], 'LineWidth', ff.lineWidth*0.75);

    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1}, 10, intStoCsensBest, 'o', 'MarkerEdgeColor', 0.5*[0.5 0.85 1], ...
    	'MarkerFaceColor', [0.5 0.85 1], 'MarkerSize', ff.markerSize*1.25, 'LineWidth', 2);
    set(theAxes{1,1}, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);

	% Finalize figure using the Publication-Ready format
    XLims = [0.1 40]; YLims = [0 1];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'eccentricity (degs)', 'Ks/Kc \times (Rs/Rc)^2');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

end
