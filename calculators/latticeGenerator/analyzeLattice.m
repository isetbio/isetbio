function analyzeLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 10; %30; 
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    %neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 32;
    
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    saveFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum));

   
    % Load results
    load(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iteration', 'maxMovements', 'terminationReason');
    fprintf('Lattice generation concluded with status: ''%s''.\n', terminationReason);
    
    visualizationParams.visualizedFOVMicrons = 200;
    
    [~, lambda] = determineMosaicWidthAndLambda(mosaicFOVDegs, neuronalType);
    
    [tabulatedDensity, tabulatedSpacing, tabulatedEcc] = ...
        generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda,  neuronalType, whichEye);
    
    setupPlotLab();
    hFig = figure(1);
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
        idx2 = find((abs(rfPositions(:,1)-350) < visualizationParams.visualizedFOVMicrons*1.4) & (absY < visualizationParams.visualizedFOVMicrons*1.4));
        idx3 = find((abs(rfPositions(:,1)-800) < visualizationParams.visualizedFOVMicrons*2.0) & (absY < visualizationParams.visualizedFOVMicrons)*2.0);
        idx4 = find((abs(rfPositions(:,1)-1400) < visualizationParams.visualizedFOVMicrons*3.0) & (absY < visualizationParams.visualizedFOVMicrons*3));
        idx5 = find((abs(rfPositions(:,1)-2400) < visualizationParams.visualizedFOVMicrons*4.0) & (absY < visualizationParams.visualizedFOVMicrons*4));
        
        plotMosaic(theAxesGrid{1,1}, rfPositions(idx,:), visualizationParams.visualizedFOVMicrons, '');
        plotMosaic(theAxesGrid{1,2}, rfPositions(idx2,:), visualizationParams.visualizedFOVMicrons*1.4, '');
        plotMosaic(theAxesGrid{1,3}, rfPositions(idx3,:), visualizationParams.visualizedFOVMicrons*2.0, '');
        plotMosaic(theAxesGrid{2,1}, rfPositions(idx4,:), visualizationParams.visualizedFOVMicrons*3.0, '');
        plotMosaic(theAxesGrid{2,2}, rfPositions(idx5,:), visualizationParams.visualizedFOVMicrons*4.0, '');
        
        triangleIndices = delaunayn(rfPositions);
        plotQuality(theAxesGrid{2,3}, rfPositions, triangleIndices, iteration, numel(maxMovements));
        drawnow;
        
        %spacingDeviations = localRFSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc);
        %plotSpacingDeviationsMap(theAxesGrid{1,2}, rfPositions(idx,:), spacingDeviations(idx), visualizationParams.visualizedFOVMicrons)
    end
end

function setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 28, ...
            'figureHeightInches', 16);
end  

