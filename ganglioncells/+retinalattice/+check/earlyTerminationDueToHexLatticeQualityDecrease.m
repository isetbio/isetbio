function [terminateNow, histogramData, widths, diffWidths, bin1Percent] = ...
    earlyTerminationDueToHexLatticeQualityDecrease(currentRFPositions, triangleIndices, widths, ...
    minHexQualityForTermination, qDistPercentile)
    
    qDist = retinalattice.compute.meshQuality(currentRFPositions, triangleIndices);
    qBins = [0.5:0.01:1.0];
    [counts,centers] = histcounts(qDist, qBins);
    bin1Percent = prctile(qDist,[qDistPercentile 3 7 15 99.8]);
    [~, idx1] = min(abs(centers-bin1Percent(2)));
    [~, idx2] = min(abs(centers-bin1Percent(3)));
    [~, idx3] = min(abs(centers-bin1Percent(4)));
    [~, idxEnd] = min(abs(centers-bin1Percent(end)));
    if (isempty(widths))
        k = 1;
    else
        k = size(widths,1)+1;
    end
    widths(k,:) = centers(idxEnd)-[centers(idx1) centers(idx2) centers(idx3)];
    if (k == 1)
        diffWidths = nan;
    else
        diffWidths = diff(widths,1)./(widths(end,:));
    end

    histogramData.x = centers;
    histogramData.y = counts;

    % Termination condition
    cond1 = bin1Percent(1) > minHexQualityForTermination;
    cond2 = (any(diffWidths(:) > 0.08)) && (~any((isnan(diffWidths(:)))));
    if (cond1 && cond2)
        terminateNow = true;
    else
        terminateNow = false;
    end
        
end
