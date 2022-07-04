function projectedCostFromOverlap = costToMaintainOverlappingInputs(obj, neighboringRGCindex, neighboringRGCconeIndices, ...
           neighboringRGCconeWeights, sourceRGCindex, sourceRGCconeIndices, sourceRGCconeWeights)

    sourceRGCConePositions = obj.inputConeMosaic.coneRFpositionsMicrons(sourceRGCconeIndices,:);
    neighboringRGCconePositions = obj.inputConeMosaic.coneRFpositionsMicrons(neighboringRGCconeIndices,:);
    
    % compute projected centroid from inputs
    sourceRGCCentroidFromInputs = RGCconnector.weightedMean(sourceRGCConePositions,sourceRGCconeWeights);
    neighborRGCCentroidFromInputs = RGCconnector.weightedMean(neighboringRGCconePositions,neighboringRGCconeWeights);

    % The more the overlap, the higher the cost
    localSpacing = 0.5*(obj.localRGCRFspacingsMicrons(neighboringRGCindex) + obj.localRGCRFspacingsMicrons(sourceRGCindex));
    d = sourceRGCCentroidFromInputs-neighborRGCCentroidFromInputs;
    separationBetweenSourceAndDestinationRGCcentroids = sqrt(sum(d.^2,2));
    
    normalizedSeparation = min([1 separationBetweenSourceAndDestinationRGCcentroids/localSpacing]);
    overlap = 1/(1+normalizedSeparation);

    % In 0-0.5 range
    projectedCostFromOverlap = (overlap-0.5);

    % Apply rfCentroidOverlapPenalty
    projectedCostFromOverlap = obj.wiringParams.rfCentroidOverlapPenaltyFactor * projectedCostFromOverlap;
end