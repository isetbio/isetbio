function plotQuality(theAxesHandle, rfPositions, triangleIndices, iteration, maxIterations)
    [minQualityValue, qValues] = computeHexLatticeQuality(rfPositions, triangleIndices);
    % Compute histogram for visualization
    qBins = 0.2:0.01:1.0;
    [histogramData.y,histogramData.x] = hist(qValues, qBins);  
    bar(theAxesHandle, histogramData.x, histogramData.y, 1, 'FaceColor', [0.5 0.5 0.5]); 
    hold(theAxesHandle, 'on');
    plot(theAxesHandle, minQualityValue*[1 1], [0 max(histogramData.y)], 'r');
    grid(theAxesHandle, 'on')
    xlabel(theAxesHandle, 'hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex');
    ylabel(theAxesHandle, '# of cells');
    %title(theAxesHandle,sprintf('iteration %d of %d', iteration, maxIterations));
end