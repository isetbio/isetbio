function visualizeCenterSurroundProperties(figNo, synthesizedRFParams, dataSet, plotlabOBJ, pdfFileName, exportsDir)
    
    % Get data
    switch (dataSet)
        case 'visual'
            rfParams = synthesizedRFParams.visual;
        case 'retinal'
             rfParams = synthesizedRFParams.retinal;
        otherwise
            error('Unknown dataSet: ''%s''.', dataSet)
    end
    eccDegs = synthesizedRFParams.eccDegs;
    
    
    hFig = figure(figNo);
    set(hFig, 'Name', sprintf('%s params', dataSet));
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.09, ...
        'widthMargin', 0.08, ...
        'heightMargin', 0.07, ...
        'bottomMargin', 0.07, ...
        'rightMargin', 0.01, ...
        'topMargin', 0.01);
    
    w = WatsonRGCModel();
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(eccDegs), ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadiiNasalRetina = 0.8 * coneRFSpacingsDegs/2;
    
    % Plot the center&surround radii
    theAxes = theAxesGrid{2,1};
    scatter(theAxes,eccDegs, rfParams.centerRadiiDegs, 'o');
    hold(theAxes, 'on');
    scatter(theAxes, eccDegs,  rfParams.surroundRadiiDegs, 'd');
    scatter(theAxes, eccDegs, coneApertureRadiiNasalRetina, 9, '.', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0]);
    %legend({'center', 'surround'});
    set(theAxes, 'XScale', 'log', 'XLim', [0.01 50], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(theAxes, 'YScale', 'log', 'YLim', [0.002 1], ...
        'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10], 'YTickLabel', {'0.001', '0.003', '0.01', '0.03',' 0.1', '0.3', '1', '3', '10'});
    xlabel(theAxes,'eccentricity (degs)');
    ylabel(theAxes,'subregion radius (degs)');
    
    % Plot the center/surround radius ratios
    theAxes = theAxesGrid{1,1};
    scatter(theAxes,eccDegs,  rfParams.centerRadiiDegs ./ rfParams.surroundRadiiDegs, 'o', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0]);
    set(theAxes, 'XLim', [0 30], 'YLim', [0 0.4], 'YTick', 0:0.05:0.4);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'center radius/surround radius');
    
    % Plot the center & surround peak sensitivities
    theAxes = theAxesGrid{2,2};
    scatter(theAxes, rfParams.centerRadiiDegs, rfParams.centerPeakSensitivities, 'o'); 
    hold(theAxes, 'on');
    scatter(theAxes,rfParams.surroundRadiiDegs, rfParams.surroundPeakSensitivities, 'd'); 
    hL = legend(theAxes, {'center', 'surround'}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 12)
    set(theAxes, 'XScale', 'log', 'XLim', [0.003 10], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30], 'XTickLabel', {'0.003', '0.01', '0.03', '0.1', '0.3', '1', '3', '10', '30'});
    set(theAxes, 'YScale', 'log', 'YLim', [0.1 40000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000]);
    xlabel(theAxes,'radius (degs)');
    ylabel(theAxes,'peak sensitivity');
    
    % Plot the integrated surround/center sensitiy ratios
    theAxes = theAxesGrid{1,2};
    integratedRatios = rfParams.surroundPeakSensitivities ./ rfParams.centerPeakSensitivities .* ...
                       (rfParams.surroundRadiiDegs ./ rfParams.centerRadiiDegs).^2;
    scatter(theAxes,synthesizedRFParams.eccDegs, integratedRatios, 'o', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0]);
    set(theAxes, 'XScale', 'log', 'XLim', [0.006 30], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30], 'YLim', [0 1]);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'integrated sensitivity (surround/center)');
    drawnow;
    
    plotlabOBJ.exportFig(hFig, 'png', sprintf('%s__%s',pdfFileName, dataSet),exportsDir);
    
end

