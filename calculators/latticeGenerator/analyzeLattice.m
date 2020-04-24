function analyzeLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 20; %30; 
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 256;
    
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    saveFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum));
   
    % Load results
    load(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iteration', 'maxMovements', 'terminationReason');
    fprintf('Lattice generation concluded with status: ''%s''.\n', terminationReason);
    
    visualizationParams.visualizedFOVMicrons = 200;
    
    setupPlotLab();
    figNo = 1;
    displayDesiredAndAchievedDensities(figNo, double(squeeze(rfPositionsHistory(end,:,:))), whichEye);
    pause
    
    
    displaySpacingDeviations = false;
    if (displaySpacingDeviations)
        [~, lambda] = determineMosaicWidthAndLambda(mosaicFOVDegs, neuronalType);
        [tabulatedDensity, tabulatedSpacing, tabulatedEcc] = ...
            generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda,  neuronalType, whichEye);
    end
    
    
    hFig = figure(2);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 3, ...
            'leftMargin', 0.04, ...
            'widthMargin', 0.07, ...
            'heightMargin', 0.07, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.06, ...
            'topMargin', 0.05);
        
    iterationsNum = size(rfPositionsHistory,1);
    for iteration = 1:iterationsNum
        % Clear axes
        for i = 1:size(theAxesGrid,1)
         for j = 1:size(theAxesGrid,2)
            cla(theAxesGrid{i,j});
         end
        end
        
        % Get the RF positions at this iteration
        rfPositions = double(squeeze(rfPositionsHistory(iteration,:,:)));
   
        % Extract the ones to visualize
        absX = abs(rfPositions(:,1));
        absY = abs(rfPositions(:,2));
        idx = find((absX < visualizationParams.visualizedFOVMicrons/2) & (absY < visualizationParams.visualizedFOVMicrons/2));
        idx2 = find((abs(rfPositions(:,1)-250) < visualizationParams.visualizedFOVMicrons*1.4) & (absY < visualizationParams.visualizedFOVMicrons*1.4));
        idx3 = find((abs(rfPositions(:,1)-600) < visualizationParams.visualizedFOVMicrons*2.0) & (absY < visualizationParams.visualizedFOVMicrons)*2.0);
        idx4 = find((abs(rfPositions(:,1)-1300) < visualizationParams.visualizedFOVMicrons*2.8) & (absY < visualizationParams.visualizedFOVMicrons*2.8));
        idx5 = find((abs(rfPositions(:,1)-2400) < visualizationParams.visualizedFOVMicrons*4.0) & (absY < visualizationParams.visualizedFOVMicrons*4));
        
        plotMosaic(theAxesGrid{1,1}, rfPositions(idx,:), visualizationParams.visualizedFOVMicrons, '');
        plotMosaic(theAxesGrid{1,2}, rfPositions(idx2,:), visualizationParams.visualizedFOVMicrons*1.4, '');
        plotMosaic(theAxesGrid{1,3}, rfPositions(idx3,:), visualizationParams.visualizedFOVMicrons*2.0, '');
        plotMosaic(theAxesGrid{2,1}, rfPositions(idx4,:), visualizationParams.visualizedFOVMicrons*3.0, '');
        plotMosaic(theAxesGrid{2,2}, rfPositions(idx5,:), visualizationParams.visualizedFOVMicrons*4.0, '');
        
        triangleIndices = delaunayn(rfPositions);
        plotQuality(theAxesGrid{2,3}, rfPositions, triangleIndices, iteration, numel(maxMovements));
        drawnow;
        
        if (displaySpacingDeviations)
            spacingDeviations = localRFSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc);
            plotSpacingDeviationsMap(theAxesGrid{1,2}, rfPositions(idx,:), spacingDeviations(idx), visualizationParams.visualizedFOVMicrons);
        end
    end
end

function displayDesiredAndAchievedDensities(figNo, positions, whichEye)
    hFig = figure(figNo);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 1, ...
            'colsNum', 3, ...
            'leftMargin', 0.04, ...
            'widthMargin', 0.05, ...
            'heightMargin', 0.01, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.05);
    
    % Instantiate 
    obj = WatsonRGCModel();

    % Compute spacings from positions
    neighborsNum = 6;
    spacingsMicrons = localRFSpacings(positions, neighborsNum);
    % Compute achieved densities from spacings
    achievedDensities = obj.densityFromSpacing(spacingsMicrons*1e-3);
    
    % Compute desired densities
    densityUnits = 'mm^2';
    [~, desiredDensities] = obj.mRGCRFSpacingAndDensityAtRetinalPositions(1e-3*positions, ...
        whichEye, 'mm', densityUnits, ...
        'adjustForISETBioConeDensity', true, 'subtype', 'ON');
    
    sampling = struct(...
        'minPos', 0, ...
        'maxPos', max(abs(positions(:))), ....
        'intervals', 100, ...
        'scale', 'linear');
    [achievedDensityMap, support] = mapFromScatteredPositions(positions, achievedDensities, sampling);
    [desiredDensityMap, support] = mapFromScatteredPositions(positions, desiredDensities, sampling);
    h = fspecial('gaussian', 5, 1);
    achievedDensityMap = conv2(achievedDensityMap, h, 'same');
    desiredDensityMap = conv2(desiredDensityMap, h, 'same');
    
    maxAchievedDensity = max(achievedDensities);
    minAchievedDensity = min(achievedDensities);
    maxDesiredDensity = max(desiredDensities);
    minDesiredDensity = min(desiredDensities);
    
    maxDensity = max([maxAchievedDensity maxDesiredDensity]);
    
    plot(theAxesGrid{1,1}, desiredDensities, achievedDensities, 'k.');
    hold(theAxesGrid{1,1}, 'on');
    plot(theAxesGrid{1,1}, [1000 maxDensity], [1000 maxDensity], 'r-');
    xlabel(theAxesGrid{1,1},'desired density (units/mm^2)')
    xlabel(theAxesGrid{1,1},'desired density (units/mm^2)')
    ylabel(theAxesGrid{1,1}, 'achieved density (units/mm^2)');
    axis(theAxesGrid{1,1}, 'square');
    set(theAxesGrid{1,1}, 'XLim', [1000 maxDensity], 'YLim', [1000 maxDensity], ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', [1000 3000 10000 30000 100000 300000], ...
        'YTick', [1000 3000 10000 30000 100000 300000], ...
        'XTickLabel', {'1k', '3k', '10k', '30k', '100k', '300k'}, ...
        'YTickLabel', {'1k', '3k', '10k', '30k', '100k', '300k'}...
        );
    grid(theAxesGrid{1,1}, 'on');
    box(theAxesGrid{1,1}, 'off');
    drawnow;
    
    xyRange = max(abs(positions(:)))*[-1 1];
    densityLevels = 1000*[75 100 125 150 200 250];
    
     
    [C,h] = contourf(theAxesGrid{1,2}, squeeze(support(:,:,1)), squeeze(support(:,:,2)), ...
        desiredDensityMap, densityLevels);
    v = densityLevels;
    clabel(C,h,v);
    
    title(theAxesGrid{1,2}, 'desired density');
    axis(theAxesGrid{1,2}, 'equal');
    set(theAxesGrid{1,2}, 'XLim', xyRange, 'YLim', xyRange);
    colormap(theAxesGrid{1,2}, brewermap(256, 'greys'));
    drawnow;
    
    [C,h] =contourf(theAxesGrid{1,3}, squeeze(support(:,:,1)), squeeze(support(:,:,2)), achievedDensityMap, densityLevels);
    v = densityLevels;
    clabel(C,h,v);
    title(theAxesGrid{1,3}, 'achieved density');
    axis(theAxesGrid{1,3}, 'equal');
    set(theAxesGrid{1,3}, 'XLim', xyRange, 'YLim', xyRange);
    colormap(theAxesGrid{1,3}, brewermap(256, 'greys'));
    
    drawnow
end


function setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'renderer', 'opengl', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 28, ...
            'figureHeightInches', 16);
end  

