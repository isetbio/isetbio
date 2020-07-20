function visualizePatchStatsDerivedFromSFcurves(patchDogModelParams, patchRGCeccentricityDegs, ...
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
    
    plotlabOBJ = setupPlotLab(0, 13, 12);
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
    
    eccRange = [0.01 30];
    radiusRange = [0.01 3];
    peakSensitivityRange = [1 100000];
    peakSensitivityRatioRange = [1e-4 1];
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
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRange, 'YTick', [1 3 10 30 100 300 1000 3000 10000 30000 100000]);
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