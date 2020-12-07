function demoCronerKaplanRGCModel
    demo1();
end

function demo1()
    ck = CronerKaplanRGCModel('generateAllFigures', true);
    
    % Generate data for the center and the surround radius
    eccDegs = logspace(log10(0.01), log10(100), 100);
    centerRadius = ck.centerRadiusFunction(ck.centerRadiusParams,eccDegs);
    centerPeakSensitivity = ck.centerPeakSensitivityFunction(ck.centerPeakSensitivityParams, centerRadius);
        
        
    p = ck.surroundRadiusParams;
    surroundRadius = ck.surroundRadiusFunction(p,eccDegs);
    surroundPeakSensitivity = ck.surroundPeakSensitivityFunction(ck.surroundPeakSensitivityParams, surroundRadius);
    
    hFig = figure(2); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
        
    
    % Plot the center radius as a function of eccentricity
    ax = theAxesGrid{1,1};
    plot(ax,eccDegs, centerRadius, '.', 'MarkerEdgeColor', [1 0 0]);
    set(ax, 'XLim', [0.01 100], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], ...
        'YLim', [0.01 10], 'YTick', [0.01 0.03 0.1 0.3 1 3 10], ...
        'XScale', 'log', 'YScale', 'log');
    ylabel(ax,'radius (degs)');
    xlabel(ax,'eccentricity (degs)');
    title(ax,'center');
    
    % Plot the surround radius as a function of eccentricity
    ax = theAxesGrid{1,2};
    plot(ax,eccDegs, surroundRadius, '.', 'MarkerEdgeColor', [0 0 1]);
    set(ax, 'XLim', [0.01 100], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], ...
        'YLim', [0.01 10], 'YTick', [0.01 0.03 0.1 0.3 1 3 10],  ...
        'XScale', 'log', 'YScale', 'log');
    xlabel(ax,'eccentricity (degs)');
    title(ax,'surround');
    
    ax = theAxesGrid{2,1};
    plot(ax,centerRadius, centerPeakSensitivity, '.', 'MarkerEdgeColor', [1 0 0]);
    set(ax, 'YLim', [0.003 3000], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000], ...
        'XLim', [0.01 10], 'XTick', [0.01 0.03 0.1 0.3 1 3 10], ...
        'XScale', 'log', 'YScale', 'log');
    xlabel(ax,'radius (degs)');
    ylabel(ax,'peak sensitivity');
  
    
    ax = theAxesGrid{2,2};
    plot(ax,surroundRadius, surroundPeakSensitivity, '.', 'MarkerEdgeColor', [0 0 1]);
    set(ax, 'YLim', [0.003 3000], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000], ...
        'XLim', [0.01 10], 'XTick', [0.01 0.03 0.1 0.3 1 3 10], ...
        'XScale', 'log', 'YScale', 'log');
    xlabel(ax,'radius (degs)');
    ylabel(ax,'peak sensitivity');

    
    
    hFig = figure(3); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
    % Plot the center radius as a function of eccentricity
    ax = theAxesGrid{1,1};
    plot(ax,eccDegs, centerRadius, '.', 'MarkerEdgeColor', [1 0 0]);
    set(ax, 'XLim', [0.01 100], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], ...
        'YLim', [0.01 10], 'YTick', [0.01 0.03 0.1 0.3 1 3 10], ...
        'XScale', 'log', 'YScale', 'log');
    ylabel(ax,'radius (degs)');
    title(ax,'center');
    
    % Plot the surround radius as a function of eccentricity
    ax = theAxesGrid{1,2};
    plot(ax,eccDegs, surroundRadius, '.', 'MarkerEdgeColor', [0 0 1]);
    set(ax, 'XLim', [0.01 100], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], ...
        'YLim', [0.01 10], 'YTick', [0.01 0.03 0.1 0.3 1 3 10],  ...
        'XScale', 'log', 'YScale', 'log');
    title(ax,'surround');
    
    % Plot the surround/center radius ratio as a function of eccentricity
    % (full ecc range)
    ax = theAxesGrid{2,1};
    plot(ax,eccDegs, centerRadius./surroundRadius, '.', 'MarkerEdgeColor', [0 0 1]);
    set(ax, 'XLim', [0 100], 'XTick', 0:10:100, ...
        'YLim', [0.01 1], 'YTick', [0.01 0.03 0.1 0.3 1], ...
        'XScale', 'linear', 'YScale', 'log');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'center radius/surround radius');
    
    % Plot the surround/center radius ratio as a function of eccentricity
    % (low eccs)
    ax = theAxesGrid{2,2};
    plot(ax,eccDegs, centerRadius./surroundRadius, '.', 'MarkerEdgeColor', [0 0 1]);
    set(ax, 'XLim', [0 5], 'XTick', 0:0.5:10, ...
        'YLim', [0 1], 'YTick', 0:0.1:1.0, ...
        'XScale', 'linear', 'YScale', 'linear');
    xlabel(ax,'eccentricity (degs)');
    
end


function demo2()
    ck = CronerKaplanRGCModel();
    
    % synthesis options - use default
    synthesisOptions = ck.synthesisOptions;
%     synthesisOptions = struct( ...
%                 'randomizeCenterRadii', ~true, ...
%                 'randomizeCenterSensitivities', ~true, ...
%                 'randomizeSurroundRadii', ~true, ...
%                 'randomizeSurroundSensitivities', true);
            
    % Simulate CronerKaplan results
    ck.simulateCronerKaplanResults('synthesisOptions', synthesisOptions, ...
        'generateVideo', true, 'generatePlots', true)
    
    

end
