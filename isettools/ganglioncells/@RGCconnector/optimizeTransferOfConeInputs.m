% Optimize how many and which of theSourceRGCinputConeIndices will
% be transfered to one of the neighboringRGCindices
function optimizeTransferOfConeInputs(obj, ...
        sourceRGC, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
        neighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights)

%    fprintf('The source RGC (%d) has %d cone inputs and %d neigbhors\n', ...
%        sourceRGC, numel(theSourceRGCinputConeIndices), numel(neighboringRGCindices));

    projectedCostsSourceRGC = zeros(numel(theSourceRGCinputConeIndices),1);
    projectedCostsNearbyRGC = zeros(numel(theSourceRGCinputConeIndices), numel(allNeighboringRGCsInputConeIndices));
    projectedCombinedCost = projectedCostsNearbyRGC;

    for iSourceConeIndex = 1:numel(theSourceRGCinputConeIndices)
        % The index of the cone to be transfered
        theConeIndex = theSourceRGCinputConeIndices(iSourceConeIndex);

        % Remove the iSourceConeIndex from theSourceRGCinputIndices
        [newSourceRGCinputConeIndices, ia] = setdiff(theSourceRGCinputConeIndices, theConeIndex);
        newSourceRGCinputConeWeights = theSourceRGCinputConeWeights(ia);

        % Compute projected cost for the source RGC
        projectedCostsSourceRGC(iSourceConeIndex) = obj.costToMaintainInputs(...
            newSourceRGCinputConeIndices, newSourceRGCinputConeWeights, ...
            obj.RGCRFspacingsMicrons(sourceRGC));

        for iNearbyRGC = 1:numel(neighboringRGCindices)
            neigboringRGC = neighboringRGCindices(iNearbyRGC);

            theNeighboringRGCinputConeIndices = allNeighboringRGCsInputConeIndices{iNearbyRGC};
            theNeighboringRGCsInputConeWeights = allNeighboringRGCsInputConeWeights{iNearbyRGC};

            % Compute the new costs for the source and nearbyRGC if the iSourceConeIndex 
            % were to be transfered from the sourceRGC to the nearbyRGC 
        
            % Add the iSourceConeIndex to theNeighboringRGCinputConeIndices
            newNeighboringRGCinputConeIndices = theNeighboringRGCinputConeIndices;
            newNeighboringRGCsInputConeWeights = theNeighboringRGCsInputConeWeights;

            newNeighboringRGCinputConeIndices(numel(newNeighboringRGCinputConeIndices)+1) = theConeIndex;
            newNeighboringRGCsInputConeWeights(numel(newNeighboringRGCsInputConeWeights)+1) = 1;

            % Compute projected cost for the nearbyRGC
            projectedCostsNearbyRGC(iSourceConeIndex, iNearbyRGC) = obj.costToMaintainInputs(...
                newNeighboringRGCinputConeIndices, newNeighboringRGCsInputConeWeights, ...
                obj.RGCRFspacingsMicrons(neigboringRGC));

            projectedCombinedCost(iSourceConeIndex, iNearbyRGC) = ...
                projectedCostsNearbyRGC(iSourceConeIndex, iNearbyRGC) + ...
                projectedCostsSourceRGC(iSourceConeIndex);
        end
    end

    % Find the (iSourceConeIndex, iNearbyRGCindex) pair that minimizes the total cost
    [~,idx] = min(projectedCombinedCost(:));
    [iSourceConeIndexOptimal, iNearbyRGCindexOptimal] = ind2sub(size(projectedCombinedCost),idx);

    % Reassign the optical cone from the sourceRGC to the optimal nearbyRGC
    obj.transferConeFromSourceRGCToDestinationRGC( ...
        theSourceRGCinputConeIndices(iSourceConeIndexOptimal), sourceRGC, neighboringRGCindices(iNearbyRGCindexOptimal));
end