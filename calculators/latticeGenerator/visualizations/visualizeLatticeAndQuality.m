function visualizeLatticeAndQuality(rfPositions, spacingDeviations, maxMovements, triangleIndices, reTriangulationIsNeeded, triangularizationTriggerEvent, iteration, iterativeParams, visualizationParams)
    hFig = figure(2); clf;
   
    if (visualizationParams.visualizeProgressOnly)
        colsNum = 1;
    else
        colsNum = 3;
    end
    theAxesGrid = plotlab.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', colsNum, ...
        'leftMargin', 0.04, ...
        'widthMargin', 0.07, ...
        'heightMargin', 0.07, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.06, ...
        'topMargin', 0.05);
    
    if (~visualizationParams.visualizeProgressOnly)
        absX = abs(rfPositions(:,1));
        absY = abs(rfPositions(:,2));
        idx = find((absX < visualizationParams.visualizedFOVMicrons/2) & (absY < visualizationParams.visualizedFOVMicrons/2));
        idx2 = find((absX < visualizationParams.visualizedFOVMicrons) & (absY < visualizationParams.visualizedFOVMicrons));

        % The central lattice on the top-left
        plotMosaic(theAxesGrid{1,2}, rfPositions(idx,:), visualizationParams.visualizedFOVMicrons, reTriangulationIsNeeded, triangularizationTriggerEvent);

        % The central density plot on the bottom-left
        %plotDensityMap(theAxesGrid{2,1}, rfPositions(idx,:), visualizedFOVMicrons);
        plotSpacingDeviationsMap(theAxesGrid{2,2}, rfPositions, spacingDeviations, visualizationParams.visualizedFOVMicrons)

        % The extended lattice on the middle
        plotMosaic(theAxesGrid{1,3}, rfPositions(idx2,:), visualizationParams.visualizedFOVMicrons*2, reTriangulationIsNeeded, triangularizationTriggerEvent);

        % The central density plot on the middle
        %plotDensityMap(theAxesGrid{2,1}, rfPositions(idx,:), visualizedFOVMicrons);
        plotSpacingDeviationsMap(theAxesGrid{2,3}, rfPositions, spacingDeviations, visualizationParams.visualizedFOVMicrons*2);
    end
    
    
    % The quality on the top-right
    plotQuality(theAxesGrid{1,1}, rfPositions, triangleIndices, iteration, iterativeParams.maxIterations);
    
    % The max movements
    plotMovements(theAxesGrid{2,1}, iteration, maxMovements, iterativeParams.dTolerance);
end

function plotSpacingDeviationsMap(theAxesHandle, rfPositions, spacingDeviations, visualizedFOVMicrons)
    
    if (~isempty(spacingDeviations))
        % Sampling vector
        sampling = struct('minPos', 1, 'maxPos', 0.5*visualizedFOVMicrons, 'intervals', 100, 'scale', 'log');
        % Generate 2D map from scattered values
        [deviationMap, mapSupport] = mapFromScatteredPositions(rfPositions, spacingDeviations, sampling);

        contourf(theAxesHandle, mapSupport(:,:,1), mapSupport(:,:,2), deviationMap, 0:0.05:1.0);
        ylabel(theAxesHandle, 'space (microns)');
        set(theAxesHandle, 'CLim', [0 1], 'ZLim', [0 1]);
        colormap(theAxesHandle,jet)
        title(theAxesHandle, sprintf('spacing deviations = %2.3f-%2.3f microns', min(deviationMap(:)), max(deviationMap(:))));
    end
    axis(theAxesHandle,'square')
end

function plotMovements(theAxesHandle, iteration, maxMovements, tolerance)
    if (~isempty(maxMovements))
        plot(theAxesHandle, 1:(iteration-1), maxMovements,'ks-'); 
        hold(theAxesHandle, 'on');
        plot(theAxesHandle, [1 iteration-1], tolerance*[1 1], 'r-');
        set(theAxesHandle, 'YScale', 'log', 'YLim', [0.0001 1], 'YTick', [0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1]);
        xlabel(theAxesHandle, 'iteration');
        ylabel(theAxesHandle, 'movement (microns)');
    end
end

function plotDensity(theAxesHandle, rfPositions, visualizedFOVMicrons)
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



function plotQuality(theAxesHandle, rfPositions, triangleIndices, iteration, maxIterations)
    [minQualityValue, qValues] = computeHexLatticeQuality(rfPositions, triangleIndices);
    % Compute histogram for visualization
    qBins = 0.2:0.01:1.0;
    [histogramData.y,histogramData.x] = hist(qValues, qBins);  
    bar(theAxesHandle, histogramData.x, histogramData.y, 1, 'FaceColor', [0.5 0.5 0.5]); 
    hold(theAxesHandle, 'on');
    plot(theAxesHandle, minQualityValue*[1 1], [0 max(histogramData.y)], 'r');
    grid(theAxesHandle, 'on')
    xlabel(theAxesHandle, 'hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex');
    title(theAxesHandle,sprintf('iteration %d of %d', iteration, maxIterations));
end


function plotMosaic(theAxesHandle, rfPositions, visualizedFOVMicrons, reTriangulationIsNeeded, triangularizationTriggerEvent)
    plot(theAxesHandle, rfPositions(:,1), rfPositions(:,2), 'k.'); 
    axis(theAxesHandle, 'equal');
    xyRange = 0.5*visualizedFOVMicrons*[-1 1];
    axis(theAxesHandle, 'square');
    set(theAxesHandle , 'XLim', xyRange, 'YLim', xyRange);
    ylabel(theAxesHandle, 'space (microns)');
    if (reTriangulationIsNeeded)
        title(theAxesHandle,sprintf('Re-triangularization in next iteration.\nReason: %s',triangularizationTriggerEvent), 'Color', 'r');
    end
end