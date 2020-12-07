function analyzeLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 20;
    
    % Type of mosaic to generate
    %neuronalType = 'cone';
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
    
    plotlabOBJ = setupPlotLab();
    figNo = 1; plotCorrespondenceAndMap = true;
    displayDesiredAndAchievedDensities(figNo, double(squeeze(rfPositionsHistory(end,:,:))), whichEye, plotCorrespondenceAndMap, plotlabOBJ);
    
    figNo = 2; plotCorrespondenceAndMap = ~true;
    displayDesiredAndAchievedDensities(figNo, double(squeeze(rfPositionsHistory(end,:,:))), whichEye, plotCorrespondenceAndMap, plotlabOBJ);
    
    figNo = 3;
    displayLatticeProgressionHistory(figNo, rfPositionsHistory, maxMovements, visualizationParams);
end

function displayLatticeProgressionHistory(figNo, rfPositionsHistory, maxMovements, visualizationParams)

    
    hFig = figure(figNo);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 4, ...
            'leftMargin', 0.02, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.03, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.05, ...
            'topMargin', 0.00);
        
    videoOBJ = VideoWriter('mosaicEvolution', 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
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
        idx4 = find((abs(rfPositions(:,1)-1300) < visualizationParams.visualizedFOVMicrons*3.5) & (absY < visualizationParams.visualizedFOVMicrons*3.5));
        idx5 = find((abs(rfPositions(:,1)-2400) < visualizationParams.visualizedFOVMicrons*5) & (absY < visualizationParams.visualizedFOVMicrons*5));
        idx6 = find((abs(rfPositions(:,1)-4000) < visualizationParams.visualizedFOVMicrons*7) & (absY < visualizationParams.visualizedFOVMicrons*7));
        
        
        plotMosaic(theAxesGrid{1,1}, rfPositions(idx,:), visualizationParams.visualizedFOVMicrons, '');
        plotMosaic(theAxesGrid{1,2}, rfPositions(idx2,:), visualizationParams.visualizedFOVMicrons*1.4, '');
        plotMosaic(theAxesGrid{1,3}, rfPositions(idx3,:), visualizationParams.visualizedFOVMicrons*2.0, '');
        plotMosaic(theAxesGrid{1,4}, rfPositions(idx4,:), visualizationParams.visualizedFOVMicrons*3.5, '');
        plotMosaic(theAxesGrid{2,1}, rfPositions(idx5,:), visualizationParams.visualizedFOVMicrons*5.0, '');
        plotMosaic(theAxesGrid{2,2}, rfPositions(idx6,:), visualizationParams.visualizedFOVMicrons*8.0, '');
        
        if (iteration<=numel(maxMovements))
            plotMovements(theAxesGrid{2,3}, iteration, maxMovements(1:iteration), []);
        end
        
        triangleIndices = delaunayn(rfPositions);
        plotQuality(theAxesGrid{2,4}, rfPositions, triangleIndices, iteration, iterationsNum);
        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
        
        displaySpacingDeviations = false;
        if (displaySpacingDeviations)
            [~, lambda] = determineMosaicWidthAndLambda(mosaicFOVDegs, neuronalType);
            [tabulatedDensity, tabulatedSpacing, tabulatedEcc] = ...
                generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda,  neuronalType, whichEye);
    
            spacingDeviations = localRFSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc);
            plotSpacingDeviationsMap(theAxesGrid{1,2}, rfPositions(idx,:), spacingDeviations(idx), visualizationParams.visualizedFOVMicrons);
        end
    end
    videoOBJ.close();
    
end

function displayDesiredAndAchievedDensities(figNo, positions, whichEye, plotCorrespondenceAndMap, plotlabOBJ)
    hFig = figure(figNo);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 1, ...
            'colsNum', 2, ...
            'leftMargin', 0.04, ...
            'widthMargin', 0.06, ...
            'heightMargin', 0.03, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.01);
    
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
    h = fspecial('gaussian', 5, 0.5);
    achievedDensityMap = conv2(achievedDensityMap, h, 'same');
    desiredDensityMap = conv2(desiredDensityMap, h, 'same');
    
    maxAchievedDensity = max(achievedDensities); minAchievedDensity = min(achievedDensities);
    maxDesiredDensity = max(desiredDensities); minDesiredDensity = min(desiredDensities);
    
    maxDensity = max([maxAchievedDensity maxDesiredDensity]);
    levels = 20;  % how many density levels
    portion = 1;  % what portion of the mosaic to display
    densityLogLevels = round(logspace(log10(300), log10(300000), levels)/100)*100;
    xyRange = max(abs(positions(:)))*[-1 1]/1.3;
    
    if (plotCorrespondenceAndMap)
        plotCorrespondence(theAxesGrid{1,1}, desiredDensities, achievedDensities, maxDensity, xyRange);
        renderContourPlot(theAxesGrid{1,2}, support, achievedDensityMap, densityLogLevels, portion*xyRange, 'achieved', true);
        drawnow;
        plotlabOBJ.exportFig(hFig, 'pdf', 'densityAndCorrespondence', pwd);
    else
        portion = 1.0;
        renderContourPlot(theAxesGrid{1,1}, support, desiredDensityMap, densityLogLevels, portion*xyRange, 'desired', true);
        renderContourPlot(theAxesGrid{1,2}, support, achievedDensityMap, densityLogLevels, portion*xyRange, 'achieved', false);
        drawnow;
        plotlabOBJ.exportFig(hFig, 'pdf', 'densityAchievedAndModel', pwd);
    end
    
    
end

function plotCorrespondence(theAxes, desiredDensities, achievedDensities, maxDensity, xyRange)
    plot(theAxes, desiredDensities, achievedDensities, 'k.');
    hold(theAxes, 'on');
    plot(theAxes, [100 maxDensity], [100 maxDensity], 'r-');
    xlabel(theAxes,'desired density (units/mm^2)');
    ylabel(theAxes, 'achieved density (units/mm^2)');
    axis(theAxes, 'equal');
    set(theAxes, 'XLim', xyRange, 'YLim', xyRange);
    set(theAxes, 'XLim', [100 maxDensity], 'YLim', [100 maxDensity], ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', [0.1 0.3 1 3 10 30 100 300]*1000, ...
        'YTick', [0.1 0.3 1 3 10 30 100 300]*1000, ...
        'XTickLabel', {'100', '300', '1k', '3k', '10k', '30k', '100k', '300k'}, ...
        'YTickLabel', {'100', '300', '1k', '3k', '10k', '30k', '100k', '300k'}...
        );
    grid(theAxes, 'on');
    box(theAxes, 'off');
end

function renderContourPlot(theAxes,support, densityMap, densityLevels, xyRange, plotTitle, showYLabel)
    [C,h] = contourf(theAxes, squeeze(support(:,:,1)), squeeze(support(:,:,2)), ...
        densityMap, densityLevels);
    v = densityLevels;
    clabel(C,h,v, 'FontSize', 16, 'LabelSpacing', 500, 'Color', 'red');
   
    title(theAxes, plotTitle);
    axis(theAxes, 'equal');
    set(theAxes, 'XLim', xyRange, 'YLim', xyRange, 'CLim', [densityLevels(1) densityLevels(end)]);
    colormap(theAxes, 0.7*brewermap(256, '*spectral')+0.3*brewermap(256, '*greys'));
    if (showYLabel)
    ylabel(theAxes, 'retinal microns');
    end
end

function plotlabOBJ = setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 28, ...
            'figureHeightInches', 14);
end  

