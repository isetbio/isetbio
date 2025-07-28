function [theSpatialCompactnessCost, theCenterConeNumerosityDifferential, theCentroidOverlapCost, theVarianceCost] = spatialCompactnessCost(...
    theTargetRFSourceRFpositions, theDonorRFSourceRFpositions, ...
    theTargetRFSourceRFconnectionWeights, theDonorRFSourceRFconnectionWeights, ...
    theTargetDestinationRFspacing, theDonorDestinationRFspacing)
	
    theTargetRFsourceRFsNum = size(theTargetRFSourceRFpositions,1);
    if (isempty(theTargetRFSourceRFconnectionWeights))
        theTargetRFSourceRFconnectionWeights = ones(theTargetRFsourceRFsNum,1);
    end
    [theTargetVariance, theTargetCentroid] = var(theTargetRFSourceRFpositions, theTargetRFSourceRFconnectionWeights, 1);
    theTargetSigma = sqrt(sum(theTargetVariance));
    
    theDonorRFsourceRFsNum = size(theDonorRFSourceRFpositions,1);
    if (isempty(theDonorRFSourceRFconnectionWeights))
        theDonorRFSourceRFconnectionWeights = ones(theDonorRFsourceRFsNum,1);
    end
    [theDonorVariance, theDonorCentroid] = var(theDonorRFSourceRFpositions, theDonorRFSourceRFconnectionWeights, 1);
    theDonorSigma = sqrt(sum(theDonorVariance));
    
	% Distance between theTargetDestinationRFindex and theDonorDestinationRFindex
    % Basically: norm(theTargetCentroid - theDonorCentroid)
    theTargetDonorDistance = sqrt(sum((theTargetCentroid - theDonorCentroid).^2,2));

    % Normalize it with respect to mean sigma of the inputs of the two destinationRFs (z-score, like)
    sigmaSum = theTargetSigma+theDonorSigma;
    theNormalizedTargetDonorSeparation = theTargetDonorDistance/sigmaSum;

	% The centroidSpacingCost is the inverse of the normalized separation between target and donor centroids
    theCentroidOverlapCost  = 1 / theNormalizedTargetDonorSeparation;

    % Add the cost of the sum of variances 
    theVarianceCost = sigmaSum / mean([theTargetDestinationRFspacing theDonorDestinationRFspacing]);

    % The centerConeNumerosityDifferential is the difference in sourceRFs between target and donor normalized by their max value
    theCenterConeNumerosityDifferential = abs(theTargetRFsourceRFsNum-theDonorRFsourceRFsNum);

    % Add the two costs
    theSpatialCompactnessCost = 0.5*(theCentroidOverlapCost + theCenterConeNumerosityDifferential);
end
