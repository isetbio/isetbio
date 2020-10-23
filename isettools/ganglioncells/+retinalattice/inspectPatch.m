function inspectPatch(patchSaveFileName)

    % Load patch generation data
    load(patchSaveFileName, 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye');
                 
    % Unpack data
    rfPositions = dataOut.rfPositions;
    rfPositionsHistory = dataOut.rfPositionsHistory;
    iterationsHistory = dataOut.iterationsHistory;
    maxMovements = dataOut.maxMovements;
    reTriangulationIterations = dataOut.reTriangulationIterations;
    terminationReason = dataOut.terminationReason;
 

    triangularizationsNum = size(reTriangulationIterations,2);
    savedIterationsNum = size(iterationsHistory,2);
    totalIterations = size(maxMovements,2);
    rfsNum = size(rfPositionsHistory,2);
    fprintf('Termination status: ''%s''.', terminationReason);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
    subplot('Position', [0.025 0.025 0.95 0.95]);
    for iter = 1:savedIterationsNum
        plot(squeeze(rfPositionsHistory(iter,:,1)), ...
             squeeze(rfPositionsHistory(iter,:,2)), 'k.');
        axis 'square'
        title(sprintf('%d', iterationsHistory(iter)));
        drawnow;
    end
    
end

