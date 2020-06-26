function plotDensityMap(theAxesHandle, rfPositions, visualizedFOVMicrons)
    % Sampling vector
    densitySampling = struct('minPos', 1/1000, 'maxPos', 0.5*visualizedFOVMicrons/1000, 'intervals', 20, 'scale', 'log');
    
    % Find distances to neighors in mm
    neighborsNum = 1;
    spacings = localRFSpacings(rfPositions/1000, neighborsNum);
    
    % Densities from medians
    densities = WatsonRGCModel.densityFromSpacing(spacings);
    
    % Generate 2D map from scattered values
    [densityMap, densityMapSupport] = mapFromScatteredPositions(rfPositions, densities, densitySampling);
    
    % plot in microns
    densityMapSupport = densityMapSupport * 1000;
    
    % Density levels to contour
    levels = [0 50 100 150 200 250 300]*1000;
    levelLabels = {'0','50k','100k','150k','200k','250k', '300k'};
    contourf(theAxesHandle, squeeze(densityMapSupport(:,:,1)), squeeze(densityMapSupport(:,:,2)),...
        densityMap, levels);
    axis(theAxesHandle, 'equal');
    xyRange = 0.5*visualizedFOVMicrons*[-1 1];
    ylabel(theAxesHandle, 'space (microns)');
    densityRange = [0 300]*1000;
    set(theAxesHandle, 'XLim', xyRange, 'YLim', xyRange, 'CLim', densityRange, 'ZLim', densityRange);
    
    % Set the color map to custom 8 color
    colormap(theAxesHandle, brewermap(1024, 'reds'));
    
    % Add colorbar without resizing the figure
%     plotlab.colorbar(theAxesHandle,'density (count/mm2)', 12, ...
%         'Location','eastoutside' ,...
%         'Ticks',levels,...
%         'TickLabels',levelLabels);
end

