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
        if (reTriangulationIsNeeded) || (strcmp(triangularizationTriggerEvent, 'final iteration'))
            mosaicTitle = sprintf('Re-triangularization in next iteration.\nReason: %s',triangularizationTriggerEvent);
        else
            mosaicTitle = '';
        end
        plotMosaic(theAxesGrid{1,2}, rfPositions(idx,:), visualizationParams.visualizedFOVMicrons, mosaicTitle);

        % The central density plot on the bottom-left
        %plotDensityMap(theAxesGrid{2,1}, rfPositions(idx,:), visualizedFOVMicrons);
        plotSpacingDeviationsMap(theAxesGrid{2,2}, rfPositions, spacingDeviations, visualizationParams.visualizedFOVMicrons)

        % The extended lattice on the middle
        plotMosaic(theAxesGrid{1,3}, rfPositions(idx2,:), visualizationParams.visualizedFOVMicrons*2, mosaicTitle);

        % The central density plot on the middle
        %plotDensityMap(theAxesGrid{2,1}, rfPositions(idx,:), visualizedFOVMicrons);
        plotSpacingDeviationsMap(theAxesGrid{2,3}, rfPositions, spacingDeviations, visualizationParams.visualizedFOVMicrons*2);
    end
    
    
    % The quality on the top-right
    plotQuality(theAxesGrid{1,1}, rfPositions, triangleIndices, iteration, iterativeParams.maxIterations);
    
    % The max movements
    plotMovements(theAxesGrid{2,1}, iteration, maxMovements, iterativeParams.dTolerance);
end


