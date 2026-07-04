% Optimize how many and which of theDestinationRFinputIndices will
% be transfered to one of the allNearbyDestinationRFindices
function optimizeTransferOfInputRFs(obj, ...
        theDestinationRFindex, theDestinationRFinputIndices, theDestinationRFinputWeights, ...
        allNearbyDestinationRFindices, allNearbyDestinationRFinputIndices, allNearbyDestinationRFinputWeights)

   %fprintf('The destinationRF (%d) has %d inputs and %d neigbhors\n', ...
   %    theDestinationRFindex, numel(theDestinationRFinputIndices), numel(allNearbyDestinationRFindices));

    projectedCostsDestinationRF = zeros(numel(theDestinationRFinputIndices),1);
    projectedCostsNearbyDestinationRF = zeros(numel(theDestinationRFinputIndices), numel(allNearbyDestinationRFinputIndices));
    projectedCombinedCost = projectedCostsNearbyDestinationRF;

    for iDestinationRFinput = 1:numel(theDestinationRFinputIndices)
        % The index of the cone to be transfered
        theDestinationRFinputIndex = theDestinationRFinputIndices(iDestinationRFinput);

        % Remove the iDestinationRFinput from theDestinationRFinputIndices
        [newDestinationRFinputIndices, ia] = setdiff(theDestinationRFinputIndices, theDestinationRFinputIndex);
        newDestinationRFinputWeights = theDestinationRFinputWeights(ia);

        % Compute projected cost for the source RGC
        theCostComponents = obj.inputMaintenanceCost(...
            newDestinationRFinputIndices, newDestinationRFinputWeights, ...
            obj.destinationRFspacingsFromCentroids(theDestinationRFindex));

        projectedCostsDestinationRF(iDestinationRFinput) = theCostComponents(1);

        for iNearbyDestinationRF = 1:numel(allNearbyDestinationRFindices)
            theNearbyDestinationRFindex = allNearbyDestinationRFindices(iNearbyDestinationRF);

            theNearbyDestinationRFinputIndices = allNearbyDestinationRFinputIndices{iNearbyDestinationRF};
            theNearbyDestinationRFinputWeights = allNearbyDestinationRFinputWeights{iNearbyDestinationRF};

            % Compute the new costs for the source and nearbyRGC if the iDestinationRFinput 
            % were to be transfered from the theDestinationRFindex to the nearbyRGC 
        
            % Add the iDestinationRFinput to theNearbyDestinationRFinputIndices
            newNearbyDestinationRFinputIndices = theNearbyDestinationRFinputIndices;
            newNearbyDestinationRFinputWeights = theNearbyDestinationRFinputWeights;

            newNearbyDestinationRFinputIndices(numel(newNearbyDestinationRFinputIndices)+1) = theDestinationRFinputIndex;
            newNearbyDestinationRFinputWeights(numel(newNearbyDestinationRFinputWeights)+1) = 1;

            % Compute projected cost for the nearbyRGC
            theCostComponents = obj.inputMaintenanceCost(...
                newNearbyDestinationRFinputIndices, newNearbyDestinationRFinputWeights, ...
                obj.destinationRFspacingsFromCentroids(theNearbyDestinationRFindex));

            projectedCostsNearbyDestinationRF(iDestinationRFinput, iNearbyDestinationRF) = ...
                theCostComponents(1);

            % Compute projected cost due to overlap
            projectedCostDueToRFoverlap = obj.overlappingDestinationRFCost(...
                              theNearbyDestinationRFindex, ...
                              newNearbyDestinationRFinputIndices, ...
                              newNearbyDestinationRFinputWeights, ...
                              theDestinationRFindex, ...
                              newDestinationRFinputIndices, ...
                              newDestinationRFinputWeights);

            projectedCombinedCost(iDestinationRFinput, iNearbyDestinationRF) = ...
                projectedCostsNearbyDestinationRF(iDestinationRFinput, iNearbyDestinationRF) + ...
                projectedCostsDestinationRF(iDestinationRFinput) + ...
                projectedCostDueToRFoverlap;
        end
    end

    % Find the (iDestinationRFinput, iNearbyDestinationRFindex) pair that minimizes the total cost
    [~,idx] = min(projectedCombinedCost(:));
    [iDestinationRFinputOptimal, iNearbyDestinationRFindexOptimal] = ind2sub(size(projectedCombinedCost),idx);

    % Reassign the optimal cone from the theDestinationRFindex to the optimal nearbyRGC
    obj.transferInputRFFromDestinationRFToNearbyDestinationRF( ...
        theDestinationRFinputIndices(iDestinationRFinputOptimal), theDestinationRFindex, allNearbyDestinationRFindices(iNearbyDestinationRFindexOptimal));
end