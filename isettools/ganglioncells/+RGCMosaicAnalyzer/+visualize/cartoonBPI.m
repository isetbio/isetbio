function cartoonBPI()
% 

	nearlyPureSurroundLMconeRatio = 0.0;
	nearlyFullMixedSurroundLMconeRatio = 0.7;
	fullMixedSurroundLMconeRatio = 1.0;

	surroundLMratios = [nearlyPureSurroundLMconeRatio fullMixedSurroundLMconeRatio]
	intSCratioMultipliers = [0.25 1 2 3 4];

	surroundLMratios = [fullMixedSurroundLMconeRatio fullMixedSurroundLMconeRatio]
	%intSCratioMultipliers = [0.25 1 2 3 4]

	centerConeColor = RGCMosaicConstructor.constants.McenterColor;
	otherConeColor = RGCMosaicConstructor.constants.LcenterColor;

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = 'demos';

	showCenterCone = true;
	showOtherCone = ~true;
	showAchromatic = ~true;
	showSTFDataPoints = true;

	for i = 1:numel(surroundLMratios)

		figNo = 1;
		[hFigScatterPlot, axScatterPlot, ff] = renderBPIaxes(figNo);

		for j = 1:numel(intSCratioMultipliers)

			figNo = 2;
			[hFigProfileAndSTF, axProfile, axSTF, ff] = renderProfileAndSTFaxes(figNo);

			[achromaticBPI(j), coneIsolatingBPI(j)] = visualizeCartoonRFprofileAndSTF(axProfile, axSTF, ff, ...
				surroundLMratios(i), intSCratioMultipliers(j), ...
				centerConeColor, otherConeColor, ...
				showCenterCone, showOtherCone, ...
				showAchromatic, showSTFDataPoints);
			set(hFigProfileAndSTF, 'Color', 'none');
			thePDFfileName = sprintf('bpiCartoonProfileAndSTF_LMratio_%2.2f_intSCratio_%2.2f.pdf', surroundLMratios(i), intSCratioMultipliers(j));
			exportBPIplot(hFigProfileAndSTF, thePDFfileName, theRawFiguresDir, pdfExportSubDir);
		end

		plot(axScatterPlot, achromaticBPI, coneIsolatingBPI, 's', 'MarkerSize', 16, ...
				'MarkerFaceColor', centerConeColor, ...
				'MarkerEdgeColor', centerConeColor*0.5);
			
		finalizeBPaxes(axScatterPlot,ff);
		thePDFfileName = sprintf('bpiCartoonScatter_LMratio_%2.2f.pdf', surroundLMratios(i));
		exportBPIplot(hFigScatterPlot, thePDFfileName, theRawFiguresDir, pdfExportSubDir);
	end
end

function [achromaticBPI, coneIsolatingBPI] = visualizeCartoonRFprofileAndSTF(axProfile, axSTF, ff, ...
	surroundLMratio, intSCratioMultiplier, centerConeColor, otherConeColor, showCenterCone, showOtherCone, showAchromatic, showSTFDataPoints)
	
	xLims = 0.3*[-1 1];
	yLims = [-0.31 1.03];

	cM = 1;
	cL = cM * surroundLMratio;
	totalLM = cL + cM;
	cL = cL/(totalLM);
	cM = cM/(totalLM);

	[cL cM]
	if (intSCratioMultiplier < 0)
		intSCratioMultiplier = 1;
	end

	eccDegs = 3.0;
	Rc = 0.02;
	intSCratio = 0.5*RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(eccDegs) * intSCratioMultiplier;
	Kc = RGCmodels.CronerKaplan.constants.centerPeakSensitivityFromCharacteristicRadiusDegsForPcells(Rc);
	Rs = Rc * RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;
	Ks = Kc * intSCratio/((Rs/Rc)^2);
	xSupport = xLims(1):0.001:xLims(2);


	% Assume center cone is L-cone
	[centerGaussianPointWeighting, centerGaussianLineWeighting, STFcenter] = gaussianRF(xSupport, Rc, Kc);
	[surroundGaussianPointWeightingLcone, surroundGaussianLineWeightingLcone, STFsurroundLcone] = gaussianRF(xSupport, Rs, Ks*cL);
	[surroundGaussianPointWeightingMcone, surroundGaussianLineWeightingMcone, STFsurroundMcone] = gaussianRF(xSupport, Rs, Ks*cM);

	achromaticLineWeightingFunction = centerGaussianLineWeighting - (surroundGaussianLineWeightingLcone+surroundGaussianLineWeightingMcone);
	achromaticSTF = STFcenter.gain - (STFsurroundLcone.gain + STFsurroundMcone.gain);
	

	Rmax = max(achromaticSTF(:));
	Ro = achromaticSTF(1);
	achromaticBPI = Ro/Rmax;
	
	coneIsolatingLineWeightingFunction = centerGaussianLineWeighting - surroundGaussianLineWeightingLcone;
	coneIsolatingSTF = STFcenter.gain - STFsurroundLcone.gain;
	Rmax = max(coneIsolatingSTF(:));
	Ro = coneIsolatingSTF(1);
	coneIsolatingBPI = Ro/Rmax;

	maxLineWeighting = max([max(achromaticLineWeightingFunction(:)) max(coneIsolatingLineWeightingFunction(:))]);
	achromaticLineWeightingFunction  = achromaticLineWeightingFunction  / maxLineWeighting;
	coneIsolatingLineWeightingFunction = coneIsolatingLineWeightingFunction / maxLineWeighting;

	hold(axProfile, 'on');
	if (showAchromatic)
		RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(axProfile, xSupport, achromaticLineWeightingFunction, 0, ...
			RGCMosaicConstructor.constants.achromaticColor, [0 0 0], 0.7, 1);
		plot(axProfile, xSupport, achromaticLineWeightingFunction, '-', 'LineWidth', 5, 'Color', RGCMosaicConstructor.constants.achromaticColor*0.5);
		plot(axProfile, xSupport, achromaticLineWeightingFunction, '-', 'LineWidth', 3, 'Color', RGCMosaicConstructor.constants.achromaticColor);
	end
	
	if (showCenterCone)
		plot(axProfile, xSupport, coneIsolatingLineWeightingFunction, '-', 'LineWidth', 5, 'Color', centerConeColor*0.5);
		plot(axProfile, xSupport, coneIsolatingLineWeightingFunction, '-', 'LineWidth', 3, 'Color', centerConeColor);
	end

	if (showOtherCone)
		otherConeIsolatingLineWeightingFunction = achromaticLineWeightingFunction - coneIsolatingLineWeightingFunction;
		plot(axProfile, xSupport, otherConeIsolatingLineWeightingFunction, '-', 'LineWidth', 5, 'Color', otherConeColor*0.5);
		plot(axProfile, xSupport, otherConeIsolatingLineWeightingFunction, '-', 'LineWidth', 3, 'Color', otherConeColor);
	end

	set(axProfile, 'XLim', xLims, 'YLim', yLims, 'Color', 'none');
	%plot(axProfile, [xLims(1) xLims(1) xLims(2) xLims(2) xLims(1)], ...
	%	     [yLims(1) yLims(2) yLims(2) yLims(1) yLims(1)], 'k-', ...
	%	     'LineWidth', 1.0);

	set(axProfile, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);
	axis(axProfile, 'square');

	% The STF
	maxSTF = max([max(coneIsolatingSTF) max(achromaticSTF)]);
	hold(axSTF, 'on');

	if (showAchromatic)
		plot(axSTF, STFcenter.sfSupport,  achromaticSTF/maxSTF, '-', 'LineWidth', 5, 'Color', RGCMosaicConstructor.constants.achromaticColor*0.3); 
		plot(axSTF, STFcenter.sfSupport,  achromaticSTF/maxSTF, '-', 'LineWidth', 3, 'Color', RGCMosaicConstructor.constants.achromaticColor);
		if (showSTFDataPoints)
			subSampledSFs = [0.1 0.3 1 2 3 5 10 20 30];
			for idx = 1:numel(subSampledSFs)
				[~, subSampledIndices(idx)] = min(abs(STFcenter.sfSupport-subSampledSFs(idx)));
			end

			plot(axSTF, STFcenter.sfSupport(subSampledIndices),  achromaticSTF(subSampledIndices)/maxSTF, 'o', 'LineWidth', 2, ...
				'MarkerSize', 20, ...
				'MarkerFaceColor', RGCMosaicConstructor.constants.achromaticColor, ...
				'MarkerEdgeColor', RGCMosaicConstructor.constants.achromaticColor*0.3); 
		end
	end
	
	if (showCenterCone)
		plot(axSTF, STFcenter.sfSupport, coneIsolatingSTF/maxSTF, '-', 'LineWidth', 5, 'Color', centerConeColor*0.3); 
		plot(axSTF, STFcenter.sfSupport, coneIsolatingSTF/maxSTF, '-', 'LineWidth', 3, 'Color', centerConeColor);

		if (showSTFDataPoints)
			subSampledSFs = [0.1 0.3 1 2 3 5 10 20 30];
			for idx = 1:numel(subSampledSFs)
				[~, subSampledIndices(idx)] = min(abs(STFcenter.sfSupport-subSampledSFs(idx)));
			end

			plot(axSTF, STFcenter.sfSupport(subSampledIndices),  coneIsolatingSTF(subSampledIndices)/maxSTF, 'o', 'LineWidth', 2, ...
				'MarkerSize', 20, ...
				'MarkerFaceColor', centerConeColor, ...
				'MarkerEdgeColor', centerConeColor*0.3); 
		end

	end
	

	if (showOtherCone)
		otherConeIsolatingSTF = abs(achromaticSTF-coneIsolatingSTF);
		plot(axSTF, STFcenter.sfSupport, otherConeIsolatingSTF/maxSTF, '-', 'LineWidth', 5, 'Color', otherConeColor*0.3); 
		plot(axSTF, STFcenter.sfSupport, otherConeIsolatingSTF/maxSTF, '-', 'LineWidth', 3, 'Color', otherConeColor);

		if (showSTFDataPoints)
			subSampledSFs = [0.1 0.3 1 2 3 5 10 20 30];
			for idx = 1:numel(subSampledSFs)
				[~, subSampledIndices(idx)] = min(abs(STFcenter.sfSupport-subSampledSFs(idx)));
			end

			plot(axSTF, STFcenter.sfSupport(subSampledIndices),  otherConeIsolatingSTF(subSampledIndices)/maxSTF, 'o', 'LineWidth', 2, ...
				'MarkerSize', 20, ...
				'MarkerFaceColor', otherConeColor, ...
				'MarkerEdgeColor', otherConeColor*0.3); 
		end

	end

	%plot(axSTF, STFcenter.sfSupport(idx), achromaticSTF(idx), 'o', 'LineWidth', 1.5, 'MarkerSize', 8, ...
	%	'MarkerFaceColor', RGCMosaicConstructor.constants.achromaticColor*0.5, 'MarkerEdgeColor', RGCMosaicConstructor.constants.achromaticColor*0.25)

	
	xLims = [0.1 100];
	yLims = [0 1.0];
	%plot(axSTF, [xLims(1) xLims(1) xLims(2) xLims(2) xLims(1)], ...
	%	     [yLims(1) yLims(2) yLims(2) yLims(1) yLims(1)], 'k-', ...
	%	     'LineWidth', 1.0);

	set(axSTF, 'XLim', xLims, 'YLim', yLims, 'XTick', [0.1 0.3 1 3 10 30 100], 'YTick', 0:0.2:1, 'YTickLabel', {}, 'XTickLabel', {});
	grid(axSTF, 'on'); box(axSTF, 'off');
	set(axSTF, 'XScale', 'log');
	%set(axSTF, 'XColor', 'none', 'YColor', 'none')
	axis(axSTF, 'square');
	PublicationReadyPlotLib.applyFormat(axSTF,ff);
	%PublicationReadyPlotLib.offsetAxes(axSTF, ff, xLims, yLims);
end

function [pointWeighting, lineWeighting, STF] = gaussianRF(x, characteristicRadius, gain)
	pointWeighting = gain * exp(-(x/characteristicRadius).^2);
	lineWeighting = pointWeighting * characteristicRadius * sqrt(pi);
	STF.sfSupport = 0.01:0.1:100;
	STF.sfSupport = [0 STF.sfSupport]
	STF.gain = gain * characteristicRadius^2 * exp(-(pi*characteristicRadius*STF.sfSupport).^2);
end

function [hFig, ax,ff] = renderBPIaxes(figNo)
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};
end

function [hFig, axProfile, axSTF, ff] = renderProfileAndSTFaxes(figNo)

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;

	axProfile = axes('Position', [0.03 0.2 0.45 0.7]);
	axSTF = axes('Position', [0.52 0.2 0.45 0.7]);
end

function finalizeBPaxes(ax,ff)
    axis (ax, 'square');
    axis(ax, 'xy');
    XLims = [0 1.0]; YLims = [0.0 1.0];
	grid(ax, 'on');
	set(ax, 'CLim', [0 1.1]);
	set(ax, 'XLim', XLims, 'XTick', 0:0.1:1, 'XTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''}, ...
		    'YLim', YLims, 'YTick', 0:0.1:1, 'YTickLabel', {'', '.1', '', '.3', '', '.5', '', '.7', '', '.9', ''});
    xlabel(ax,'achromatic gratings');
    ylabel(ax,'cone-isolating gratings');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);
end

function exportBPIplot(hFig, pdfFileName, theRawFiguresDir, pdfExportSubDir)

	thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end

