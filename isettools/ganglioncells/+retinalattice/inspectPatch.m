function inspectPatch(fovDegs, neuronType, whichEye)

    % Configure algorithm params
    params = retinalattice.configure(fovDegs, neuronType, whichEye);
    
    % Load patch generation data
    load(fullfile(params.latticeGalleryDir, params.patchSaveFileName), 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye');
                 

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
    
    figure(2); clf;
    plot(1:size(maxMovements,2), maxMovements, 'ko-');
    set(gca, 'YScale', 'log');
    
        
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
    subplot('Position', [0.025 0.025 0.95 0.95]);
    
    xyLimsMicrons = 1000*[-0.5 0.5]*RGCmodels.Watson.convert.rhoDegsToMMs(fovDegs);
    for iter = 1:savedIterationsNum
    
        theRFpositions = double(squeeze(rfPositionsHistory(iter,:,:)));
        triangleIndices = retinalattice.compute.triangularization(theRFpositions, ...
            params.domainFunction, params.radius, params.borderTolerance);
   
        qDist = retinalattice.compute.meshQuality(theRFpositions, triangleIndices);
        bin1Percent = prctile(qDist,params.qDistPercentile);
        subplot(1,2,1);
        plot(theRFpositions(:,1), theRFpositions(:,2), 'k.');
        set(gca, 'XLim', xyLimsMicrons, 'YLim', xyLimsMicrons);
        axis 'square';
        title(sprintf('%d', iterationsHistory(iter)));
        
        subplot(1,2,2);
        qBins = [0.5:0.01:1.0];
        [counts,centers] = histcounts(qDist, qBins);
        bar(centers(2:end), counts); hold on;
        plot(bin1Percent(1)*[1 1], [0 max(counts)], 'r-', 'LineWidth', 1.5);
        set(gca, 'XLim', [qBins(1) qBins(end)], 'YLim', [0 max(counts)]);
        title(sprintf('Iter %d', iterationsHistory(iter)));
        axis 'square';
        hold off
        drawnow;
        
       
    end
    
end

