function visualizeRFparamsForConnectedPatch(figNo, figName, ...
    RGCconeInputInfo, ...
    eccRadiusDegs, ...
    centerCharacteristicRadii, surroundCharacteristicRadii, ...
    centerPeakSensitivities, surroundPeakSensitivities, ...
    pdfFigureFileName, figExportsDir, plotlabOBJ)

    % Ranges and ticks
    eccRange = [0.3 30];
    eccTicks = [0.1 0.3 1 3 10 30 100];
    
    radiusRange = [0.003 3];
    radiusTicks = [0.003 0.01 0.03 0.1 0.3 1 3 10];
    
    peakSensitivityRange = [0.1 3000];
    peakSensitivityTicks = [0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000];
    
    peakSensitivityRatioRange = [1e-4 1];
    peakSensitivityRatioTicks = [0.0001 0.001 0.01  0.1  1];
    
    % Compute cone aperture radius in eccRange
    w = WatsonRGCModel();
    eccCones = logspace(log10(0.3), log10(30), 100);
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(eccCones, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadiiNasalRetina = WatsonRGCModel.coneApertureToDiameterRatio * 0.5 * coneRFSpacingsDegs;
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(eccCones, ...
            'superior meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadiiSuperiorRetina = WatsonRGCModel.coneApertureToDiameterRatio * 0.5 * coneRFSpacingsDegs;
    
    ck = CronerKaplanRGCModel('generateAllFigures', false, 'instantiatePlotLab', false);
    
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
    
    % The cell eccentricities
    ax = theAxesGrid{1,1};hold(ax, 'on');
    scatter(ax, eccRadiusDegs, eccRadiusDegs*0, ...
        144, 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none');
    scatter(ax, ck.recordedEccentricitiesData(:,1), ck.recordedEccentricitiesData(:,2), 121);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'eccentricity (degs)');
    set(ax, 'XLim', [-30 30] , 'XTick', -30:10:30);
    set(ax,  'YLim',  [-30 30] , 'YTick', -30:10:30);
    axis(ax, 'square');
    
    % center characteristic radii as a function of eccentricity
    ax = theAxesGrid{1,2}; hold(ax, 'on');
    scatter(ax, eccRadiusDegs, centerCharacteristicRadii, ...
        144, 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none'); 
    scatter(ax, ck.centerData('size').eccDegs, ck.centerData('size').radiusDegs, ...
        121); 
    line(ax, eccCones, coneApertureRadiiNasalRetina, 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    line(ax, eccCones, coneApertureRadiiSuperiorRetina, 'LineStyle', ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'center radius (degs)');
    set(ax, 'XScale', 'log', 'XLim', eccRange , 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', radiusRange, 'YTick',radiusTicks);
    axis(ax, 'square');
    

    % surround characteristic radii as a function of eccentricity
    ax = theAxesGrid{1,3}; hold(ax, 'on');
    scatter(ax, eccRadiusDegs, surroundCharacteristicRadii, ...
        144, 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none'); 
    scatter(ax, ck.surroundData('size').eccDegs, ck.surroundData('size').radiusDegs, 121);
    
    line(ax, eccCones, coneApertureRadiiNasalRetina, 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    line(ax, eccCones, coneApertureRadiiSuperiorRetina, 'LineStyle', ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'surround radius (degs)');
    set(ax, 'XScale', 'log', 'XLim', eccRange , 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', radiusRange, 'YTick',radiusTicks);
    axis(ax, 'square');
    
    
    
     % ratio of characteristic radii (S:C)
    ax = theAxesGrid{2,1}; hold(ax, 'on');
    characteristicRadiiRatios = surroundCharacteristicRadii./centerCharacteristicRadii;
    scatter(ax,eccRadiusDegs, characteristicRadiiRatios,...
        144, 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none'); 
    scatter(ax, ck.surroundCenterRatiosData('radius').eccDegs, ...
        ck.surroundCenterRatiosData('radius').ratios, 121);

    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'characteristic radius (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'linear', 'YLim', [0 25], 'YTick', 0:5:25);
    axis(ax, 'square');
    
    % ratio of peak sensitivities (S:C)
    ax = theAxesGrid{2,2}; hold(ax, 'on');
    peakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;
    scatter(ax,eccRadiusDegs, peakSensitivityRatios, ...
        144, 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none'); 
    scatter(ax, ck.surroundCenterRatiosData('peak sensitivity').eccDegs, ...
        ck.surroundCenterRatiosData('peak sensitivity').ratios, 121); 
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'peak sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'log', 'YLim', peakSensitivityRatioRange, 'YTick', peakSensitivityRatioTicks);
    axis(ax, 'square');
    
   
    
    % ratio of integrated sensitivities (S:C)
    ax = theAxesGrid{2,3}; hold(ax, 'on');
    integratedSensitivityRatios = (surroundCharacteristicRadii./centerCharacteristicRadii).^2 .* (surroundPeakSensitivities./centerPeakSensitivities);
    scatter(ax, eccRadiusDegs, integratedSensitivityRatios, ...
        144, 'MarkerEdgeColor', 'none');
    scatter(ax, ck.surroundCenterRatiosData('integrated sensitivity').eccDegs, ...
        ck.surroundCenterRatiosData('integrated sensitivity').ratios, 121);
    
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'integrated sensitivity (surround/center)');
    set(ax, 'XScale', 'log', 'XLim', eccRange, 'XTick', eccTicks);
    set(ax, 'YScale', 'linear', 'YLim', [0 1], 'YTick', 0:0.2:1);
    axis(ax, 'square');
    
    if (~isempty(pdfFigureFileName))
         plotlabOBJ.exportFig(hFig, 'pdf', pdfFigureFileName, figExportsDir);
    end
    
    exportDataFile = sprintf('%s.mat', fullfile(figExportsDir, pdfFigureFileName));
    save(exportDataFile, 'eccRadiusDegs', 'LMconeBalance', ...
        'centerCharacteristicRadii', 'surroundCharacteristicRadii', ...
        'centerPeakSensitivities', 'surroundPeakSensitivities');
end

