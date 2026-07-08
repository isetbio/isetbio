function theTotalPoolingCosts = totalPoolingCosts(obj)

    destinationRFsNum = size(obj.connectivityMatrix,2);
    theTotalConeInputUniformityCosts = zeros(1, destinationRFsNum);
    theTotalCenterConeNumerosityDifferentials = zeros(1,destinationRFsNum);
    theTotalCentroidOverlapCosts = zeros(1,destinationRFsNum);
    theTotalVarianceCosts = zeros(1,destinationRFsNum);
    theTotalSpatialUniformityCosts = zeros(1, destinationRFsNum);

    parfor theTargetDestinationRF = 1:destinationRFsNum
        theTargetDestinationRFsourceRFindices = find(squeeze(obj.connectivityMatrix(:,theTargetDestinationRF))>0.001);
        theTargetDestinationRFsourceRFpositions = obj.sourceLattice.RFpositionsMicrons(theTargetDestinationRFsourceRFindices,:);
        inputConeTypes = obj.sourceLattice.metaData.coneTypes(theTargetDestinationRFsourceRFindices);
        theTargetDestinationRFsourceWeights = full(obj.connectivityMatrix(theTargetDestinationRFsourceRFindices,theTargetDestinationRF));
        theTotalConeInputUniformityCosts(theTargetDestinationRF) = coneToMidgetRGCConnector.spectralUniformityCost(inputConeTypes, theTargetDestinationRFsourceWeights);

        % Find the indices of the neigboring destinationRFs. These come sorted in increasing distance from theTargetDestinationRF
        [theNearbyDestinationRFindicesToTheTargetDestinationRF, theDistancesToTheTargetDestinationRF] = ...
            obj.indicesOfNeighboringDestinationRFs(...
                theTargetDestinationRF, ...
                'maxNeighborsNum', 6);

        nearbyDestinationRFsNum = numel(theNearbyDestinationRFindicesToTheTargetDestinationRF);
        theConeNumerosityDifferentialForThisPairOfNearbyDestinationRFs = zeros(1,nearbyDestinationRFsNum);
        theCentroidOverlapCostForThisPairOfNearbyDestinationRFs = zeros(1,nearbyDestinationRFsNum);
        theVarianceCostForThisPairOfNearbyDestinationRFs = zeros(1,nearbyDestinationRFsNum);

        for iNearbyDestinationRF = 1:nearbyDestinationRFsNum
            theNearbyDestinationRFindex = theNearbyDestinationRFindicesToTheTargetDestinationRF(iNearbyDestinationRF);
            theNearbyDestinationRFinputSourceRFindices = find(squeeze(obj.connectivityMatrix(:,theNearbyDestinationRFindex))>0.001);
            theNearbyDestinationRFsourceRFpositions = obj.sourceLattice.RFpositionsMicrons(theNearbyDestinationRFinputSourceRFindices,:);
            theNearbyDestinationRFsourceWeights = full(obj.connectivityMatrix(theNearbyDestinationRFinputSourceRFindices,theNearbyDestinationRFindex));

            [~, theConeNumerosityDifferentialForThisPairOfNearbyDestinationRFs(iNearbyDestinationRF), ...
            theCentroidOverlapCostForThisPairOfNearbyDestinationRFs(iNearbyDestinationRF), ...
            theVarianceCostForThisPairOfNearbyDestinationRFs(iNearbyDestinationRF)] = ...
                coneToMidgetRGCConnector.spatialCompactnessCost(...
                    theTargetDestinationRFsourceRFpositions, theNearbyDestinationRFsourceRFpositions, ...
                    theTargetDestinationRFsourceWeights, theNearbyDestinationRFsourceWeights, ...
                    obj.destinationRFspacingsFromCentroids(theTargetDestinationRF), obj.destinationRFspacingsFromCentroids(theNearbyDestinationRFindex));

        end % for iNearbyDestinationRF

        theTotalCenterConeNumerosityDifferentials(theTargetDestinationRF) = ...
            sum(theConeNumerosityDifferentialForThisPairOfNearbyDestinationRFs) / nearbyDestinationRFsNum;

        theTotalCentroidOverlapCosts(theTargetDestinationRF) = ...
            sum(theCentroidOverlapCostForThisPairOfNearbyDestinationRFs)/ nearbyDestinationRFsNum;

        theTotalVarianceCosts(theTargetDestinationRF) = ...
            sum(theVarianceCostForThisPairOfNearbyDestinationRFs)/ nearbyDestinationRFsNum;

    end  %for theTargetDestinationRF

    theTotalConeInputUniformityCost = sum(theTotalConeInputUniformityCosts) / destinationRFsNum;
    theTotalCenterConeNumerosityDifferential = sum(theTotalCenterConeNumerosityDifferentials)  / destinationRFsNum;
    theTotalCentroidOverlapCost = sum(theTotalCentroidOverlapCosts) / destinationRFsNum;
    theTotalVarianceCost = sum(theTotalVarianceCosts)/destinationRFsNum;

    % The total cost: weighted sum of theTotalSpatialUniformityCost and theTotalConeInputUniformityCost
    theTotalPoolingCosts(1,1) = ...
        obj.wiringParams.spatialChromaticUniformityTradeoff * (theTotalCenterConeNumerosityDifferential + theTotalCentroidOverlapCost) + ...
        (1-obj.wiringParams.spatialChromaticUniformityTradeoff) * theTotalConeInputUniformityCost;
    theTotalPoolingCosts(1,2) = theTotalCenterConeNumerosityDifferential;
    theTotalPoolingCosts(1,3) = theTotalConeInputUniformityCost;
    theTotalPoolingCosts(1,4) = theTotalCentroidOverlapCost;
    theTotalPoolingCosts(1,5) = theTotalVarianceCost;

end