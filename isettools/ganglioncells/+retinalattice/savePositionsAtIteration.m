% Visualizes the progression of the mosaic (quality and max movement) and
% saves the lattice positions at the desired iteration
function savePositionsAtIteration(fovDegs, neuronType, whichEye)

    % Configure algorithm params
    p = retinalattice.configure(fovDegs, neuronType, whichEye);
    mosaicProgressFileName = fullfile(p.latticeGalleryDir, p.patchSaveFileName);
    mosaicFinalFileName = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
    % Load patch generation data
    load(mosaicProgressFileName, 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye');

    allIterations = 1:size(dataOut.maxMovements,2);
    figure(); clf;
    subplot(1,2,1);
    plot(allIterations, dataOut.maxMovements, 'ko-');
    ylabel('max movement');
    set(gca, 'YScale', 'log');
    subplot(1,2,2);
    plot(dataOut.iterationsHistory, dataOut.qualityHistory, 'ko-');
    ylabel('quality');
    drawnow;
    
    % Query user about iteration no.
    iter = input('Iteration to save: ');
    ii = find((dataOut.iterationsHistory-iter)<=0);
    targetIteration = ii(end);

    % Extract rfPositions from history
    rfPositions = squeeze(dataOut.rfPositionsHistory(targetIteration,:,:));
    
    % Export
    save(mosaicFinalFileName, ...
        'rfPositions', 'params', 'fovDegs', 'neuronType', 'whichEye', '-v7.3');
    fprintf('Positions from iteration %d saved in %s\n', dataOut.iterationsHistory(targetIteration), ...
        mosaicFinalFileName);
end

