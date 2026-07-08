function [theSpatialCompactnessCosts, theSpectralUniformityCosts, ...
          atLeastOneConeInputWithGreaterThanMinWeightIncluded, theInputConeNumerosityDifferentials, ...
          theCentroidOverlapCosts, theSpatialVarianceCosts] = ...
        rfCenterSpatioChromaticCosts(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.parse(varargin{:});
    minConeWeightIncluded = p.Results.minConeWeightIncluded;
    
    % Initialize
    theSpectralUniformityCosts = zeros(1, obj.rgcsNum);
    theSpatialCompactnessCosts = zeros(1, obj.rgcsNum);
    theCentroidOverlapCosts  = zeros(1, obj.rgcsNum);
    theSpatialVarianceCosts  = zeros(1, obj.rgcsNum);
    theInputConeNumerosityDifferentials = zeros(1, obj.rgcsNum);

    atLeastOneConeInputWithGreaterThanMinWeightIncluded = false(1, obj.rgcsNum);

    centerConnectivityMatrix = obj.rgcRFcenterConeConnectivityMatrix;
    coneRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons;
    coneTypes = obj.inputConeMosaic.coneTypes;

    

    parfor iRGC = 1:obj.rgcsNum
        connectivityVector = squeeze(centerConnectivityMatrix(:, iRGC));
        inputConeIndices = find(connectivityVector > minConeWeightIncluded);
        
        theTargetRGCinputConeTypes = coneTypes(inputConeIndices);
        theTargetRGCinputConeWeights = full(connectivityVector(inputConeIndices));
        theTargetRGCinputConeNumerosity = numel(theTargetRGCinputConeWeights);
        theTargetRGCinputConePositions = coneRFpositionsMicrons(inputConeIndices,:);

        if (numel(inputConeIndices) > 0)
            atLeastOneConeInputWithGreaterThanMinWeightIncluded(iRGC) = true;
        end

        theSpectralUniformityCosts(iRGC) = coneToMidgetRGCConnector.spectralUniformityCost(...
            theTargetRGCinputConeTypes, theTargetRGCinputConeWeights);
        
        % Find the indices of the neigboring destinationRFs. These come sorted in increasing distance from theTargetDestinationRF
        [theNearbyRGCindicesToTheTargetRGC, theDistancesToTheTargetRGC] = ...
            obj.indicesOfNearbyRGCs(iRGC, ...
                'maxNeighborsNum', 6);

        nearbyDestinationRFsNum = numel(theNearbyRGCindicesToTheTargetRGC);
        theNearbySpatialCompactnessCosts = zeros(1,nearbyDestinationRFsNum);
        theNearbyCentroidOverlapCosts  = zeros(1,nearbyDestinationRFsNum);
        theNearbySpatialVarianceCosts = zeros(1,nearbyDestinationRFsNum);

        theNearbyRGCinputConeNumerosityDifferentials = zeros(1,nearbyDestinationRFsNum);
        for iNearbyRGC = 1:nearbyDestinationRFsNum
            theNearbyRGCindex = theNearbyRGCindicesToTheTargetRGC(iNearbyRGC);
            connectivityVector = squeeze(centerConnectivityMatrix(:, theNearbyRGCindex));
            theNearbyRGCinputConeIndices = find(connectivityVector >  minConeWeightIncluded);

            theNearbyRGCinputConeWeights = full(connectivityVector(theNearbyRGCinputConeIndices));
            theNearbyRGCinputConeNumerosityDifferentials(iNearbyRGC) = abs(numel(theNearbyRGCinputConeWeights)-theTargetRGCinputConeNumerosity);
            theNearbyRGCinputConePositions = coneRFpositionsMicrons(theNearbyRGCinputConeIndices,:);

            theTargetDestinationRFspacing = [];
            theNearbyDestinationRFspacing = [];

            switch (class(obj))
                case 'MosaicConnector'
                    theTargetDestinationRFspacing = obj.destinationRFspacingsFromCentroids(iRGC);
                    theNearbyDestinationRFspacing = obj.destinationRFspacingsFromCentroids(theNearbyRGCindex);
                case  'mRGCMosaic'
                    theTargetDestinationRFspacing = obj.rgcRFspacingsMicrons(iRGC);
                    theNearbyDestinationRFspacing = obj.rgcRFspacingsMicrons(theNearbyRGCindex);
                otherwise
                    error('Unknown class')
            end

            [theNearbySpatialCompactnessCosts(iNearbyRGC), ~, theNearbyCentroidOverlapCosts(iNearbyRGC), theNearbySpatialVarianceCosts(iNearbyRGC)] = ...
                coneToMidgetRGCConnector.spatialCompactnessCost(theTargetRGCinputConePositions, theNearbyRGCinputConePositions, ...
                    theTargetRGCinputConeWeights, theNearbyRGCinputConeWeights, ...
                    theTargetDestinationRFspacing, theNearbyDestinationRFspacing);
        end % for iNearbyRGC

        theSpatialCompactnessCosts(iRGC) = mean(theNearbySpatialCompactnessCosts);
        theCentroidOverlapCosts(iRGC) = mean(theNearbyCentroidOverlapCosts);
        theSpatialVarianceCosts(iRGC) = mean(theNearbySpatialVarianceCosts);
        theInputConeNumerosityDifferentials(iRGC) = max(theNearbyRGCinputConeNumerosityDifferentials);
    end % iRGC

end
