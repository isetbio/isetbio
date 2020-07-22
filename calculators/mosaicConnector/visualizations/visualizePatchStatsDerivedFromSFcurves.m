function visualizePatchStatsDerivedFromSFcurves(patchDogModelParams, synthesizedRFParams, patchRGCeccentricityDegs, ...
        LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir)
    
    % Preallocate memory
    rgcsNum = numel(patchDogModelParams);
    centerRadii = zeros(1, rgcsNum);
    surroundRadii = zeros(1, rgcsNum);
    centerPeakSensitivities = zeros(1, rgcsNum);
    surroundPeakSensitivities = zeros(1, rgcsNum);
    
    % Extract params from the fitted model
    for iRGC = 1:rgcsNum
        p = patchDogModelParams{iRGC};
        centerRadii(iRGC) = p.rC;
        surroundRadii(iRGC) = p.rS;
        centerPeakSensitivities(iRGC) = p.kC;
        surroundPeakSensitivities(iRGC) = p.kS;
    end
    
    eccRange = [0.01 30];
    radiusRange = [0.01 3];
    peakSensitivityRange = [0.1 3000];
    peakSensitivityRatioRange = [1e-4 1];
    
    plotlabOBJ = setupPlotLab(0, 13, 12);
    hFig = figure(554); clf
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.03, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.08, ...
        'bottomMargin', 0.07, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    ax = theAxesGrid{1,1};
    scatter(ax, synthesizedRFParams.visual.centerRadiiDegs, centerRadii); hold(ax, 'on');
    scatter(ax, synthesizedRFParams.visual.surroundRadiiDegs, surroundRadii);
    line(ax, radiusRange, radiusRange, 'Color', 'k', 'LineWidth', 1.5);
    set(ax, 'XScale', 'log', 'XLim', radiusRange , 'XTick', [0.01 0.03 0.1 0.3 1 3]);
    set(ax, 'YScale', 'log', 'YLim', radiusRange, 'YTick', [0.01 0.03 0.1 0.3 1 3]);
    axis(ax, 'square');
    xlabel(ax,'synthesized radius (degs)');
    ylabel(ax,'fitted model radius (degs)');
    
    ax = theAxesGrid{1,2};
    scatter(ax, synthesizedRFParams.visual.centerPeakSensitivities, centerPeakSensitivities); hold(ax, 'on');
    scatter(ax, synthesizedRFParams.visual.surroundPeakSensitivities, surroundPeakSensitivities);
    line(ax, peakSensitivityRange, peakSensitivityRange, 'Color', 'k', 'LineWidth', 1.5);
    set(ax, 'XScale', 'log', 'XLim', peakSensitivityRange, 'XTick', [1 3 10 30 100 300 1000 3000 10000 30000 100000]);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRange, 'YTick', [1 3 10 30 100 300 1000 3000 10000 30000 100000]);
    axis(ax, 'square');
    xlabel(ax,'synthesized peak sensitivity');
    ylabel(ax,'fitted peak sensitivity');
   
    ax = theAxesGrid{2,1};
    synthesizedIntegratedSensitivityRatios = (synthesizedRFParams.visual.surroundPeakSensitivities .* synthesizedRFParams.visual.surroundRadiiDegs.^2) ./ ...
                                             (synthesizedRFParams.visual.centerPeakSensitivities .* synthesizedRFParams.visual.centerRadiiDegs.^2);
    sensitivityRatios =  (surroundPeakSensitivities .* surroundRadii.^2) ./ ...
                         (centerPeakSensitivities .* centerRadii.^2);        
    scatter(ax, synthesizedIntegratedSensitivityRatios, sensitivityRatios, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5]); hold(ax, 'on');
    line(ax, [0 1] , [0 1], 'Color', 'k', 'LineWidth', 1.5);
    set(ax, 'XScale', 'linear', 'XLim', [0 1], 'XTick', 0:0.2:1);
    set(ax, 'YScale', 'linear', 'YLim', [0 1], 'YTick', 0:0.2:1);
    axis(ax, 'square');
    xlabel(ax, 'synthesized integr. sensitivity ratio (surround/center)');
    ylabel(ax, 'fitted integr. sensitivity ratio (surround/center)');
    
    ax = theAxesGrid{2,2};
    synthesizedRadiusRatios = synthesizedRFParams.visual.centerRadiiDegs ./synthesizedRFParams.visual.surroundRadiiDegs;
    radiusRatios = centerRadii ./surroundRadii;
    scatter(ax, synthesizedRadiusRatios, radiusRatios, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5]); hold(ax, 'on');
    line(ax, [0 1] , [0 1], 'Color', 'k', 'LineWidth', 1.5);
    set(ax, 'XScale', 'linear', 'XLim', [0 1], 'XTick', 0:0.2:1);
    set(ax, 'YScale', 'linear', 'YLim', [0 1], 'YTick', 0:0.2:1);
    axis(ax, 'square');
    xlabel(ax, 'synthesized radius ratio (center/surround)');
    ylabel(ax, 'fitted radius ratio (center/surround)');
    
    if (exportFig)
         pdfFileName = sprintf('PatchStatsVsSynthesizedStats_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s', LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
         plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, figExportsDir);
    end
    pause
    
    
    hFig = figure(555); clf
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.03, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.08, ...
        'bottomMargin', 0.07, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    
    ax = theAxesGrid{1,1};
    scatter(ax, patchRGCeccentricityDegs, centerRadii, 'MarkerEdgeColor', 'none'); hold(ax, 'on');
    scatter(ax, patchRGCeccentricityDegs, surroundRadii, 'MarkerEdgeColor', 'none');
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'radius (degs)');
    
    set(ax, 'XScale', 'log', 'XLim', eccRange , 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'log', 'YLim', radiusRange, 'YTick', [0.01 0.03 0.1 0.3 1 3]);
    axis(ax, 'square');
    
    ax = theAxesGrid{1,2};
    scatter(ax, centerRadii,  centerPeakSensitivities, 'MarkerEdgeColor', 'none' ); hold(ax, 'on');
    scatter(ax, surroundRadii, surroundPeakSensitivities, 'MarkerEdgeColor', 'none' );
    ylabel(ax, 'peak sensitivity');
    xlabel(ax, 'radius (degs)');
    set(ax, 'XScale', 'log', 'XLim', radiusRange, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRange, 'YTick', [0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000]);
    axis(ax, 'square');

    ax = theAxesGrid{2,1};
    peakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;
    scatter(ax,patchRGCeccentricityDegs, peakSensitivityRatios , 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5]);
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'peak sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRatioRange, 'YTick', [0.0001 0.001 0.01  0.1  1]);
    axis(ax, 'square');
    
    ax = theAxesGrid{2,2};
    integratedSensitivityRatios = (surroundRadii/centerRadii).^2 .* (surroundPeakSensitivities./centerPeakSensitivities);
    scatter(ax, patchRGCeccentricityDegs, integratedSensitivityRatios, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5]);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'integrated sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'linear', 'YLim', [0 1], 'YTick', [0:0.2:1]);
    axis(ax, 'square');

    if (exportFig)
         pdfFileName = sprintf('PatchStats_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s', LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
         plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, figExportsDir);
    end
    setupPlotLab(-1);
end


function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 