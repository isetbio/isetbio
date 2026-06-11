function theCostComponents = overlappingDestinationRFCost(obj, ...
            theRGCindex, ...
            theRGCinputIndices, theRGCinputWeights, ...
            theNearbyRGCindex, ...
            theNearbyRGCinputIndices, theNearbyRGCinputWeights ...
            )

    sourceRGCConePositions = obj.sourceLattice.RFpositionsMicrons(theRGCinputIndices,:);
    neighboringRGCconePositions = obj.sourceLattice.RFpositionsMicrons(theNearbyRGCinputIndices,:);
    
    % compute projected centroid from inputs
    sourceRGCCentroidFromInputs = MosaicConnector.weightedMean(sourceRGCConePositions,theRGCinputWeights);
    neighborRGCCentroidFromInputs = MosaicConnector.weightedMean(neighboringRGCconePositions,theNearbyRGCinputWeights);

    % The more the overlap, the higher the cost
    localSpacing = 0.5*(obj.destinationRFspacingsFromCentroids(theRGCindex) + ...
                        obj.destinationRFspacingsFromCentroids(theNearbyRGCindex));
    d = sourceRGCCentroidFromInputs-neighborRGCCentroidFromInputs;
    separationBetweenSourceAndDestinationRGCcentroids = sqrt(sum(d.^2,2));
    
    normalizedSeparation = min([1 separationBetweenSourceAndDestinationRGCcentroids/localSpacing]);
    overlap = 1/(1+normalizedSeparation);

    % In 0-0.5 range
    projectedCostFromOverlap = (overlap-0.5);

    % Apply rfCentroidOverlapPenalty
    theCostComponents = obj.wiringParams.rfCentroidOverlapPenaltyFactor * projectedCostFromOverlap;
end

