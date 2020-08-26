function visualizeRFparamsForConnectedPatch(figNo, figName, ...
    RGCconeInputInfo, ...
    eccRadiusDegs, ...
    centerCharacteristicRadii, surroundCharacteristicRadii, ...
    centerPeakSensitivities, surroundPeakSensitivities, ...
    pdfFigureFileName, figExportsDir, plotlabOBJ)

    % Ranges and ticks
    eccRange = [0.01 30];
    eccTicks = [0.01 0.03 0.1 0.3 1 3 10 30];
    
    radiusRange = [0.001 1];
    radiusTicks = [0.01 0.03 0.1 0.3 1 3];
    
    peakSensitivityRange = [0.1 3000];
    peakSensitivityTicks = [0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000];
    
    peakSensitivityRatioRange = [1e-4 1];
    peakSensitivityRatioTicks = [0.0001 0.001 0.01  0.1  1];
    
    % Compute cone aperture radius in eccRange
    w = WatsonRGCModel();
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(eccRadiusDegs), ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadiiNasalRetina = WatsonRGCModel.coneApertureToDiameterRatio * 0.5 * coneRFSpacingsDegs;
    
    % Determine RF center cone dominance
    LMconeBalance = determineConeDominance(RGCconeInputInfo);
    
    hFig = figure(figNo); clf
    set(hFig, 'Name', figName);
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 3, ...
        'leftMargin', 0.03, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.08, ...
        'bottomMargin', 0.07, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    if (isempty(LMconeBalance))
        markerColorsCenter = [1 0 0];
        markerColorsSurround = [0 0 1];
    else
        markerColorsCenter = zeros(numel(centerCharacteristicRadii),3);
        markerColorsCenter(:,1) = LMconeBalance.center;
        markerColorsCenter(:,2) = 1-LMconeBalance.center;
        markerColorsCenter(:,3) = 0.3;
        markerColorsSurround = zeros(numel(centerCharacteristicRadii),3);
        markerColorsSurround(:,1) = LMconeBalance.surround;
        markerColorsSurround(:,2) = 1-LMconeBalance.surround;
        markerColorsSurround(:,3) = 0.3;
    end
    
    % center and surround characteristic radii as a function of eccentricity
    ax = theAxesGrid{1,1};
    scatter(ax, eccRadiusDegs, centerCharacteristicRadii, [], markerColorsCenter, 'o', 'MarkerEdgeColor', 'none'); hold(ax, 'on');
    scatter(ax, eccRadiusDegs, surroundCharacteristicRadii, [], markerColorsSurround, 'd', 'MarkerEdgeColor', 'none');
    line(ax, eccRadiusDegs, coneApertureRadiiNasalRetina, 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'radius (degs)');
    set(ax, 'XScale', 'log', 'XLim', eccRange , 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', radiusRange, 'YTick',radiusTicks);
    axis(ax, 'square');
    
    % center and surround peak sensitivities as a function of eccentricity
    ax = theAxesGrid{1,2};
    scatter(ax, eccRadiusDegs, centerPeakSensitivities, [], markerColorsCenter, 'o', 'MarkerEdgeColor', 'none'); hold(ax, 'on');
    scatter(ax, eccRadiusDegs, surroundPeakSensitivities, [], markerColorsSurround, 'd', 'MarkerEdgeColor', 'none');
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'peak sensitivity');
    set(ax, 'XScale', 'log', 'XLim', eccRange , 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRange, 'YTick', peakSensitivityTicks);
    axis(ax, 'square');
    
    % center and surround peak sensitivities as a function of respective characteristic radii
    ax = theAxesGrid{1,3};
    scatter(ax, centerCharacteristicRadii,  centerPeakSensitivities, [], markerColorsCenter, 'o', 'MarkerEdgeColor', 'none' ); hold(ax, 'on');
    scatter(ax, surroundCharacteristicRadii, surroundPeakSensitivities, [], markerColorsCenter, 'd', 'MarkerEdgeColor', 'none' );
    ylabel(ax, 'peak sensitivity');
    xlabel(ax, 'radius (degs)');
    set(ax, 'XScale', 'log', 'XLim', radiusRange, 'XTick', radiusTicks);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRange, 'YTick', peakSensitivityTicks);
    axis(ax, 'square');
    
    % ratio of peak sensitivities (S:C)
    ax = theAxesGrid{2,1};
    peakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;
    scatter(ax,eccRadiusDegs, peakSensitivityRatios , [], markerColorsCenter, 'MarkerEdgeColor', 'none');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'peak sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRatioRange, 'YTick', peakSensitivityRatioTicks);
    axis(ax, 'square');
    
    % ratio of characteristic radii (S:C)
    ax = theAxesGrid{2,2};
    characteristicRadiiRatios = surroundCharacteristicRadii./centerCharacteristicRadii;
    scatter(ax,eccRadiusDegs, characteristicRadiiRatios, [], markerColorsCenter, 'MarkerEdgeColor', 'none');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'characteristic radius (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'linear', 'YLim', [0 10], 'YTick', 0:2:10);
    axis(ax, 'square');
    
    % ratio of integrated sensitivities (S:C)
    ax = theAxesGrid{2,3};
    integratedSensitivityRatios = (surroundCharacteristicRadii./centerCharacteristicRadii).^2 .* (surroundPeakSensitivities./centerPeakSensitivities);
    scatter(ax, eccRadiusDegs, integratedSensitivityRatios, [], markerColorsCenter, 'MarkerEdgeColor', 'none');
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'integrated sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'linear', 'YLim', [0 1], 'YTick', 0:0.2:1);
    axis(ax, 'square');
    
    if (~isempty(pdfFigureFileName))
         plotlabOBJ.exportFig(hFig, 'pdf', pdfFigureFileName, figExportsDir);
    end
    
end

