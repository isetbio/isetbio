% Phase 7: Visualize weights of cone inputs to mRGC  RF centers and surrounds
function runPhase7(runParams)
    % Weights datafile
    weightsFile = fullfile(runParams.outputDir, sprintf('%s_inPatchAt_%2.1f_%2.1fdegs_WithSize_%2.2f_%2.2f.mat', ...
        runParams.inputFile, runParams.patchEccDegs(1), runParams.patchEccDegs(2), runParams.patchSizeDegs(1), runParams.patchSizeDegs(2)));
    
    load(weightsFile, ...
        'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios', ...
        'midgetRGCconnectionMatrixCenter', 'midgetRGCconnectionMatrixSurround', ...
        'synthesizedRFParams', 'patchEccDegs', 'patchSizeDegs');
    
    % Set up plotlab
    plotlabOBJ = setupPlotLab();
    outputFile = sprintf('%s_Summary',runParams.outputFile);
    
    % Visualize summary of visual rf properties
    %visualizeCenterSurroundProperties(1, synthesizedRFParams, 'visual', plotlabOBJ, outputFile,runParams.exportsDir);
    
    % Visualize summary of retinal rf properties
    %visualizeCenterSurroundProperties(2, synthesizedRFParams, 'retinal', plotlabOBJ, outputFile,runParams.exportsDir);
    
    % Visualize # of cones and cone types per RF center as a function of eccentricity
    visualizeConeInputRatioToSubregions(3, midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
        synthesizedRFParams.rgcIndices, synthesizedRFParams.eccDegs, coneTypes, ...
        plotlabOBJ, outputFile, runParams.exportsDir);
    pause
    
    % Visualize examples of retinal 2D RFs (video)
    outputFile = sprintf('%s_RFexamples',runParams.outputFile);
    visualizeSubregions(4,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
        synthesizedRFParams.rgcIndices,  synthesizedRFParams.eccDegs, ...
        synthesizedRFParams.centerPositionMicrons, synthesizedRFParams.retinal.centerRadiiDegs, synthesizedRFParams.retinal.surroundRadiiDegs,...
        conePositionsMicrons, coneSpacingsMicrons,  coneTypes, ...
        plotlabOBJ, outputFile, runParams.exportsDir);
    
end


function visualizeConeInputRatioToSubregions(figNo, midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
        rgcIndices, rgcEccDegs, coneTypes, plotlabOBJ, pdfFileName,exportsDir)
    
    hFig = figure(figNo); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.08, ...
        'widthMargin', 0.08, ...
        'heightMargin', 0.1, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.01, ...
        'topMargin', 0.01);
    
    binEdges = 0:0.5:30;
    inputs1Count = zeros(1, numel(binEdges));
    inputs2Count = zeros(1, numel(binEdges));
    inputs3Count = zeros(1, numel(binEdges));
    inputs4Count = zeros(1, numel(binEdges));
    inputs2CountSameCone = zeros(1, numel(binEdges));
    inputs3CountSameCone = zeros(1, numel(binEdges));
    inputs4CountSameCone = zeros(1, numel(binEdges));
    
    for iRGC  = 1:numel(rgcIndices)   
        
        % Get the rgcIndex (index in the full RGC mosaic, not in the
        % patch we are currently analyzing)
        rgcIndex = rgcIndices(iRGC);
        
        % center weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixCenter(:, rgcIndex)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        centerWeights = squeeze(full(midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, rgcIndex)));
        
        iecc = round(rgcEccDegs(iRGC)/0.5)+1;
        switch (numel(coneIndicesConnectedToCenter))
            case 1
                inputs1Count(iecc) = inputs1Count(iecc) + 1;
            case 2
                inputs2Count(iecc) = inputs2Count(iecc) + 1;
                if (coneTypes(coneIndicesConnectedToCenter(1)) == coneTypes(coneIndicesConnectedToCenter(2)))
                    inputs2CountSameCone(iecc) = inputs2CountSameCone(iecc) +1;
                end
            case 3
                inputs3Count(iecc) = inputs3Count(iecc) + 1;
                if (coneTypes(coneIndicesConnectedToCenter(1)) == coneTypes(coneIndicesConnectedToCenter(2))) && ...
                   (coneTypes(coneIndicesConnectedToCenter(1)) == coneTypes(coneIndicesConnectedToCenter(3)))
                    inputs3CountSameCone(iecc) = inputs3CountSameCone(iecc) +1;
                end
            otherwise
                inputs4Count(iecc) = inputs4Count(iecc) + 1;
                allSame = true;
                for k = 1:numel(coneIndicesConnectedToCenter)
                    if (coneTypes(coneIndicesConnectedToCenter(k))~= coneTypes(coneIndicesConnectedToCenter(1)))
                        allSame = false;
                    end
                end
                if (allSame)
                    inputs4CountSameCone(iecc) = inputs4CountSameCone(iecc) +1;
                end
                
        end
        
        
    end % iRGC
    
    binEdges = binEdges + 0.25;
    theAxes = theAxesGrid{1,1};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes,  binEdges, inputs1Count, 1, 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-0.25, inputs1Count, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', [0 25], 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', [0 3000]);
    ylabel(theAxes,'number of  rgcs');
    title(theAxes,'single cone center')
    
    theAxes = theAxesGrid{1,2};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges,  inputs2Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-0.25, inputs2CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', [0 25], 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', [0 150]);
    title(theAxes,'2 cone center')
    
    theAxes = theAxesGrid{2,1};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges,  inputs3Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-0.25, inputs3CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', [0 25], 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', [0 150]);
    xlabel(theAxes,'eccentricity (degs)');
    ylabel(theAxes,'number of rgcs');
    title(theAxes,'3 cone center')
    
    theAxes = theAxesGrid{2,2};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges, inputs4Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-0.25, inputs4CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', [0 25], 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', [0 150]);
    xlabel(theAxes,'eccentricity (degs)');
    title(theAxes,'4+ cone center')
    
    plotlabOBJ.exportFig(hFig, 'png', sprintf('%s__numberOfConesInCenter',pdfFileName),exportsDir);
    
    
end

function visualizeSubregions(figNo,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, rgcIndices, rgcEccDegs, ...
    rgcCenterPositionMicrons, retinalCenterRadiiDegs, retinalSurroundRadiiDegs, conePositionsMicrons, coneSpacingsMicrons,  coneTypes, plotlabOBJ, pdfFileName,exportsDir)

    maxGainDivisor = 500;

    hFig = figure(figNo); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.08, ...
        'widthMargin', 0.08, ...
        'heightMargin', 0.02, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.01, ...
        'topMargin', 0.01);
    
   
    videoOBJ = VideoWriter(fullfile(exportsDir, pdfFileName), 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();
        
    
    previousEccDegs = [];
    previousCenterRadii = [];
    previousSurroundRadii = [];
    previousIntegratedSensitivityRatios = [];
    
    randomizedK = randperm(numel(rgcIndices));
    for k = 1:numel(randomizedK)   
        
        % Get randomized index
        iRGC = randomizedK(k);
        
        % Get the rgcIndex (index in the full RGC mosaic, not in the
        % patch we are currently analyzing)
        rgcIndex = rgcIndices(iRGC);
        
        % Get position of RGC in microns
        rgcPositionMicrons = rgcCenterPositionMicrons(iRGC,:);
         
        % center weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixCenter(:, rgcIndex)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        centerWeights = squeeze(full(midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, rgcIndex)));
        
        % surround weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixSurround(:, rgcIndex)));
        coneIndicesConnectedToSurround = find(connectivityVector>0);
        surroundWeights = squeeze(full(midgetRGCconnectionMatrixSurround(coneIndicesConnectedToSurround, rgcIndex)));
        
        
        % Get radius of RGC center in microns
        rgcCenterRadiusDegs = retinalCenterRadiiDegs(iRGC);
        rgcCenterRadiusMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(rgcCenterRadiusDegs, rgcEccDegs(iRGC));
        
        % Compute visualized spatial support
        rgcSurroundRadiusDegs = retinalSurroundRadiiDegs(iRGC);
        rgcSurroundRadiusMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(rgcSurroundRadiusDegs, rgcEccDegs(iRGC));
        spatialSupportRangeMicrons = rgcSurroundRadiusMicrons*2; 
        spatialSupportRangeMicrons = round(0.2*300);
        
        % Plot center weights
        % Only visualize weights up to a fraction of the mean(center). After that saturate to max
        maxWeightVisualized = max(centerWeights);
        
        visualizeWeigtedConeInputsToRGCSubregion(theAxesGrid{1,1}, coneTypes(coneIndicesConnectedToCenter), ...
            conePositionsMicrons(coneIndicesConnectedToCenter,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToCenter), ...
            centerWeights, maxWeightVisualized, ...
            rgcPositionMicrons, rgcCenterRadiusMicrons, spatialSupportRangeMicrons, 1);
        
        % Plot surround weights
        % Only visualize weights up to a fraction of the mean(center). After that saturate to max
        maxWeightVisualized = max(surroundWeights);
        visualizeWeigtedConeInputsToRGCSubregion(theAxesGrid{1,2}, coneTypes(coneIndicesConnectedToSurround), ...
            conePositionsMicrons(coneIndicesConnectedToSurround,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToSurround), ...
            surroundWeights,  maxWeightVisualized, ...
            rgcPositionMicrons, rgcSurroundRadiusMicrons, spatialSupportRangeMicrons, -1);
        
        % Plot the center&surround radii
        theAxes = theAxesGrid{2,1};
        cla(theAxes);
        hold(theAxes, 'on');
        if (~isempty(previousEccDegs))
             scatter(theAxes, previousEccDegs, previousCenterRadii, 'o');
             scatter(theAxes, previousEccDegs, previousSurroundRadii, 'd');
        end
        
        scatter(theAxes, rgcEccDegs(iRGC), rgcCenterRadiusDegs, 225,'o', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 0 0]);
        scatter(theAxes, rgcEccDegs(iRGC),  rgcSurroundRadiusDegs, 225,'d', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [0 0 1]);
        
        set(theAxes, 'XScale', 'log', 'XLim', [0.006 30], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30]);
        set(theAxes, 'YScale', 'log', 'YLim', [0.002 1], ...
            'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10], 'YTickLabel', {'0.001', '0.003', '0.01', '0.03',' 0.1', '0.3', '1', '3', '10'});
        xlabel(theAxes,'eccentricity (degs)');
        ylabel(theAxes,'subregion radius (degs)');
    
        % Plot integrated sensitivity ratio
        theAxes = theAxesGrid{2,2};
        cla(theAxes);
        hold(theAxes, 'on');
        if (~isempty(previousEccDegs))
            scatter(theAxes, previousEccDegs, previousIntegratedSensitivityRatios, 'o', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0]);
        end
        scatter(theAxes, rgcEccDegs(iRGC), sum(surroundWeights)/sum(centerWeights), 225, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [0 0 0]);
        
        set(theAxes, 'XScale', 'log', 'XLim', [0.006 30], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30], 'YLim', [0 1]);
        xlabel(theAxes, 'eccentricity (degs)');
        ylabel(theAxes, 'integrated sensitivity (surround/center)');

        previousEccDegs = cat(2, previousEccDegs, rgcEccDegs(iRGC));
        previousCenterRadii = cat(2, previousCenterRadii, rgcCenterRadiusDegs);
        previousSurroundRadii = cat(2, previousSurroundRadii, rgcSurroundRadiusDegs);
        previousIntegratedSensitivityRatios = cat(2, previousIntegratedSensitivityRatios, sum(surroundWeights)/sum(centerWeights));
        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
         
        %plotlabOBJ.exportFig(hFig, 'png', pdfFileName,exportsDir);
    end
    videoOBJ.close();
end


function plotlabOBJ = setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'lineMarkerSize', 6, ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 14, ...
            'figureHeightInches', 14);
end 

