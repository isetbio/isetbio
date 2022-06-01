function projectedCostFromOverlap = costToMaintainOverlappingInputs(obj, neighboringRGCindex, neighboringRGCconeIndices, ...
           neighboringRGCconeWeights, sourceRGCindex, sourceRGCconeIndices, sourceRGCconeWeights)

    sourceRGCConePositions = obj.inputConeMosaic.coneRFpositionsMicrons(sourceRGCconeIndices,:);
    neighboringRGCconePositions = obj.inputConeMosaic.coneRFpositionsMicrons(neighboringRGCconeIndices,:);
    
    % compute projected centroid from inputs
    [~, sourceRGCCentroidFromInputs] = var(sourceRGCConePositions,sourceRGCconeWeights,1);
    [~, neighborRGCCentroidFromInputs] = var(neighboringRGCconePositions,neighboringRGCconeWeights,1);

    % The more the overlap, the higher the cost
    localSpacing = 0.5*(obj.localRGCRFspacingsMicrons(neighboringRGCindex) + obj.localRGCRFspacingsMicrons(sourceRGCindex));
    d = sourceRGCCentroidFromInputs-neighborRGCCentroidFromInputs;
    separationBetweenSourceAndDestinationRGCcentroids = sqrt(sum(d.^2,2));
    
    normalizedSeparation = min([1 separationBetweenSourceAndDestinationRGCcentroids/localSpacing]);
    overlap = 1/(1+normalizedSeparation);
    projectedCostFromOverlap = 2*(overlap-0.5);
end