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
    
    allIterations = 1:1:size(maxMovements,2);
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1600 700], 'Color', [1 1 1]);
    
   
    subplot(1,2,1);
    plot(allIterations, maxMovements, 'ko-');
    hold on;
    set(gca, 'YScale', 'log');
    for iter = 1:savedIterationsNum
        plot(iterationsHistory(iter)*[1 1], maxMovements(iterationsHistory(iter)), 'ro', 'MarkerFaceColor', [1 0.5 0.5]);
        drawnow
        theRFpositions = double(squeeze(rfPositionsHistory(iter,:,:)));
        triangleIndices = retinalattice.compute.triangularization(theRFpositions, ...
            params.domainFunction, params.radius, params.borderTolerance);
   
        qDist{iter} = retinalattice.compute.meshQuality(theRFpositions, triangleIndices);
        bin1Percent(iter) = prctile(qDist{iter},params.qDistPercentile);
    end
    xlabel('iteration no');
    ylabel('movement');
    set(gca, 'FontSize', 14);
    
    
    subplot(1,2,2);
    plot(iterationsHistory, bin1Percent, 'rs-', 'LineWidth', 1.5);
    set(gca, 'YLim', [0.5 1.0]);
    xlabel('iteration no');
    ylabel('q value');
    set(gca, 'FontSize', 14);
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1600 700], 'Color', [1 1 1]);
    
    xyLimsMicrons = 1000*[-0.5 0.5]*RGCmodels.Watson.convert.rhoDegsToMMs(fovDegs);
    for iter = 1:savedIterationsNum
    
        theRFpositions = double(squeeze(rfPositionsHistory(iter,:,:)));
        theQDist = qDist{iter};
        subplot(1,2,1);
        plot(theRFpositions(:,1), theRFpositions(:,2), 'k.');
        set(gca, 'XLim', xyLimsMicrons, 'YLim', xyLimsMicrons,  'FontSize', 14);
        axis 'square';
        title(sprintf('%d', iterationsHistory(iter)));
        
        subplot(1,2,2);
        qBins = [0.5:0.01:1.0];
        [counts,centers] = histcounts(theQDist, qBins);
        bar(centers(2:end), counts); hold on;
        plot(bin1Percent(iter)*[1 1], [0 max(counts)], 'r-', 'LineWidth', 1.5);
        set(gca, 'XLim', [qBins(1) qBins(end)], 'YLim', [0 max(counts)], 'FontSize', 14);
        xlabel('q value');
        title(sprintf('Iter %d', iterationsHistory(iter)));
        axis 'square';
        hold off
        drawnow;
        
       
    end
    
end

