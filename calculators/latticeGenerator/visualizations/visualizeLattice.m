function visualizeLattice(rfPositions)

    % Compute density in count/mm2
    densitySampling = struct('minPos', 1/1000, 'maxPos', max(abs(rfPositions(:)))/1000, 'intervals', 20, 'scale', 'log');
    [densityMap, densityMapSupport] = densityMapFromPositions(rfPositions/1000, densitySampling);
    % plot in microns
    densityMapSupport = densityMapSupport * 1000;
    
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 18, ...
            'figureHeightInches', 9);
        
    hFig = figure(1);
   
    theAxesGrid = plotlab.axesGrid(hFig, ...
        'rowsNum', 1, ...
        'colsNum', 2, ...
        'leftMargin', 0.04, ...
        'widthMargin', 0.07, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.06, ...
        'topMargin', 0.05);
    
    % The lattice on the left
    theAxesHandle = theAxesGrid{1,1};
    scatter(theAxesHandle, rfPositions(:,1), rfPositions(:,2), 25);
    axis(theAxesHandle, 'equal');
    xyRange = max(abs(rfPositions(:)))*[-1 1];
    set(theAxesHandle , 'XLim', xyRange, 'YLim', xyRange);
    plotlab.offsetAxes(theAxesHandle, 'offsetPercent', 0.015);
    
    % The density map on the right
    theAxesHandle = theAxesGrid{1,2};
    levels = [0 50 100 150 200 250 300]*1000;
    levelLabels = {'0','50k','100k','150k','200k','250k', '300k'};
    contourf(squeeze(densityMapSupport(:,:,1)), squeeze(densityMapSupport(:,:,2)),...
        densityMap, levels);
    axis(theAxesHandle, 'equal');
    xRange = max(max(abs(squeeze(densityMapSupport(:,:,1)))))*[-1 1];
    yRange = max(max(abs(squeeze(densityMapSupport(:,:,2)))))*[-1 1];
    densityRange = [0 300]*1000;
    set(theAxesHandle, 'XLim', xRange, 'YLim', yRange, 'CLim', densityRange, 'ZLim', densityRange);
    
    % Set the color map to custom 8 color
    colormap(theAxesHandle, brewermap(1024, 'reds'));
    
    % Add colorbar without resizing the figure
    plotlab.colorbar(theAxesHandle,'density (count/mm2)', 12, ...
        'Location','eastoutside' ,...
        'Ticks',levels,...
        'TickLabels',levelLabels);
    
end
    
