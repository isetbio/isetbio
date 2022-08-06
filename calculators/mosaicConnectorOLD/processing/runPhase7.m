% Phase 7: Visualize weights of cone inputs to mRGC  RF centers and surrounds
function runPhase7(runParams)
    % Weights datafile
    
    % Assemble cone weights data file name
    qNum = numel(runParams.deconvolutionOpticsParams.quadrantsToAverage);
    quadrantsToAverage = '';
    for qIndex = 1:qNum
        quadrantsToAverage = sprintf('_%s', quadrantsToAverage, runParams.deconvolutionOpticsParams.quadrantsToAverage{qIndex});
    end
    
    dataFile = sprintf('%s_inPatchAt_%2.1f_%2.1fdegs_WithSize_%2.2f_%2.2f_ForSubject_%d_AndQuadrants%s.mat', ...
                        runParams.inputFile, runParams.patchEccDegs(1), runParams.patchEccDegs(2), runParams.patchSizeDegs(1), runParams.patchSizeDegs(2), ...
                        runParams.deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage, ...       // Deconvolution model: which subject)
                        quadrantsToAverage...                                            // Deconvolution model: which quadrant to use/average
                );    
    
    load(fullfile(runParams.outputDir, dataFile), ...
        'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios', ...
        'midgetRGCconnectionMatrixCenter', 'midgetRGCconnectionMatrixSurround', ...
        'synthesizedRFParams', 'patchEccDegs', 'patchSizeDegs', ...
        'PolansSubjectIDsAveraged', 'quadrantsAveraged');
    
    % Set up plotlab
    plotlabOBJ = setupPlotLab();
    outputFile = sprintf('%s_Summary',runParams.outputFile);

    % Visualize summary of visual rf properties
    %visualizeCenterSurroundProperties(1, synthesizedRFParams, 'visual', plotlabOBJ, outputFile,runParams.exportsDir);
    
    % Visualize summary of retinal rf properties
    %visualizeCenterSurroundProperties(2, synthesizedRFParams, 'retinal', plotlabOBJ, outputFile,runParams.exportsDir);
    
    % Visualize # of cones and cone types per RF center as a function of eccentricity
%     visualizeConeInputRatioToSubregions(3, midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
%         synthesizedRFParams.eccDegs, coneTypes, ...
%         plotlabOBJ, outputFile, runParams.exportsDir);
%     pause
    
    % Visualize examples of retinal 2D RFs (video)
    outputFile = sprintf('%s_RFexamples',runParams.outputFile);
    visualizeSubregions(4,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
         synthesizedRFParams.eccDegs, ...
        synthesizedRFParams.centerPositionMicrons, synthesizedRFParams.retinal.centerRadiiDegs, synthesizedRFParams.retinal.surroundRadiiDegs,...
        conePositionsMicrons, coneSpacingsMicrons,  coneTypes, ...
        plotlabOBJ, outputFile, runParams.exportsDir);
    
end


function visualizeConeInputRatioToSubregions(figNo, midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
        rgcEccDegs, coneTypes, plotlabOBJ, pdfFileName,exportsDir)
    
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
    
    eccBinWidthDegs = 0.25;
    binEdges = 0:eccBinWidthDegs:30;
    binEdges = binEdges + eccBinWidthDegs/2;
    
    inputs1Count = zeros(1, numel(binEdges));
    inputs2Count = zeros(1, numel(binEdges));
    inputs3Count = zeros(1, numel(binEdges));
    inputs4Count = zeros(1, numel(binEdges));
    inputs2CountSameCone = zeros(1, numel(binEdges));
    inputs3CountSameCone = zeros(1, numel(binEdges));
    inputs4CountSameCone = zeros(1, numel(binEdges));
    
    rgcsNum = size(midgetRGCconnectionMatrixCenter,2);
    for rgcIndex  = 1:rgcsNum   
        
        % center weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixCenter(:, rgcIndex)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        centerWeights = squeeze(full(midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, rgcIndex)));
        
        iecc = round(rgcEccDegs(rgcIndex)/eccBinWidthDegs)+1;
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
        
        
    end % rgcIndex
    
    eccRange = [0 28];
    cellsRange = [0 12000];
    cellsRange2 = [0 6000];
    theAxes = theAxesGrid{1,1};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes,  binEdges, inputs1Count, 1, 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-eccBinWidthDegs/4, inputs1Count, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', eccRange, 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', cellsRange);
    ylabel(theAxes,'number of  rgcs');
    title(theAxes,'single cone center')
    
    theAxes = theAxesGrid{1,2};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges,  inputs2Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-eccBinWidthDegs/4, inputs2CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', eccRange, 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', cellsRange2);
    title(theAxes,'2 cone center')
    
    theAxes = theAxesGrid{2,1};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges,  inputs3Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-eccBinWidthDegs/4, inputs3CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', eccRange, 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', cellsRange2);
    xlabel(theAxes,'eccentricity (degs)');
    ylabel(theAxes,'number of rgcs');
    title(theAxes,'3 cone center')
    
    theAxes = theAxesGrid{2,2};
    cla(theAxes);
    hold(theAxes, 'on');
    bar(theAxes, binEdges, inputs4Count, 1,  'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none'); 
    stairs(theAxes, binEdges-eccBinWidthDegs/4, inputs4CountSameCone, 'k-', 'LineWidth', 1.5);
    set(theAxes, 'XScale', 'linear', 'XLim', eccRange, 'XTick', 0:5:25);
    set(theAxes, 'YScale', 'linear', 'YLim', cellsRange2);
    xlabel(theAxes,'eccentricity (degs)');
    title(theAxes,'4+ cone center')
    
    plotlabOBJ.exportFig(hFig, 'png', sprintf('%s__numberOfConesInCenter',pdfFileName),exportsDir);
    
    
end


