function convergenceIsAchieved = convergenceAchieved(netTotalCostSequence)

    trackingValues = 3;
    currentPass = numel(netTotalCostSequence);
    if (currentPass >= 2*trackingValues)
        initialIndices = (-trackingValues): -1 : (-2*trackingValues+1);
        finalIndices = 0 : -1 : (-trackingValues+1);
        initialValue = sum(netTotalCostSequence(numel(netTotalCostSequence)+initialIndices));
        finalValue = sum(netTotalCostSequence(numel(netTotalCostSequence)+finalIndices));
        epsilon = 100*abs(finalValue-initialValue)/finalValue;
    else
        epsilon = inf;
    end
    
    
    if (epsilon < 1)
        convergenceIsAchieved = true;
    else
        convergenceIsAchieved = false;
    end
    
end