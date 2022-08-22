function beneficialSwapWasFound = optimizeSwappingOfInputRFs(obj,...
        theDestinationRFindex, theDestinationRFinputIndices, theDestinationRFinputWeights, ...
        allNearbyDestinationRFindices, allNearbyDestinationRFinputIndices, allNearbyDestinationRFinputWeights)

    % Compute cost for the destinationRF to maintain its current inputs
    theDestinationRFCostComponents = obj.inputMaintenanceCost(...
        theDestinationRFinputIndices, theDestinationRFinputWeights, ...
        obj.destinationRFspacingsFromCentroids(theDestinationRFindex));

    
    % Only swap up to a max of the inputs -1
    theDestinationRFinputsNum = numel(theDestinationRFinputIndices);
    maxDestinationInputsNumSwapped = max([1 ceil(theDestinationRFinputsNum/2)]);

    % Cost for all combinations of nearby destination RFs and # of inputs to be swapped
    optimalProjectedCosts = inf(numel(allNearbyDestinationRFindices), 100);
    optimalDestinationRFinputCombination = cell(numel(allNearbyDestinationRFindices), 100);
    optimalNearbyDestinationRFinputCombination = cell(numel(allNearbyDestinationRFindices), 100);

    % Costs before swapping
    costsBeforeSwapping = zeros(1, numel(allNearbyDestinationRFindices));
 
    for iNearbyDestinationRF = 1:numel(allNearbyDestinationRFindices)

        % Compute combined source + neighbor RGC cost before swapping
        theNearbyDestinationRFindex = allNearbyDestinationRFindices(iNearbyDestinationRF);

        theNearbyDestinationRFinputIndices = allNearbyDestinationRFinputIndices{iNearbyDestinationRF};
        theNearbyDestinationRFinputWeights = allNearbyDestinationRFinputWeights{iNearbyDestinationRF};

        theNearbyRFCostComponents = obj.inputMaintenanceCost(theNearbyDestinationRFinputIndices, ...
                                         theNearbyDestinationRFinputWeights, ...
                                         obj.destinationRFspacingsFromCentroids(theNearbyDestinationRFindex));
        
        costsBeforeSwapping(iNearbyDestinationRF) = ...
                theDestinationRFCostComponents(1) + ...
                theNearbyRFCostComponents(1) + ...
                obj.overlappingDestinationRFCost(...
                                         theNearbyDestinationRFindex, ...
                                         theNearbyDestinationRFinputIndices, ...
                                         theNearbyDestinationRFinputWeights, ...
                                         theDestinationRFindex, ...
                                         theDestinationRFinputIndices, ...
                                         theDestinationRFinputWeights);

        % Find max # of input RFs to be swapped
        nearbyDestinationRFinputsNum = numel(theNearbyDestinationRFinputIndices);
        maxNearbyDestinationInputsNumSwapped = max([1 ceil(nearbyDestinationRFinputsNum/2)]);

        maxInputsNumSwapped = min([...
            maxDestinationInputsNumSwapped ...
            maxNearbyDestinationInputsNumSwapped ...
            obj.wiringParams.maxNumberOfConesToSwap]);

        for inputsNumSwapped = 1:maxInputsNumSwapped
            % Compute all inputsNumSwapped permutations for the source RGC cone inputs
            destinationRFinputCombinations = nchoosek(1:theDestinationRFinputsNum, inputsNumSwapped);

            % Compute all inputsNumSwapped permutations for the neighboring RGC cone inputs
            nearbyDestinationRFinputCombinations = nchoosek(1:nearbyDestinationRFinputsNum, inputsNumSwapped);

            % Allocate memory to store all the projected costs
            destinationRFinputCombinationsNum = size(destinationRFinputCombinations,1);
            nearbyDestinationRFinputCombinationsNum = size(nearbyDestinationRFinputCombinations,1);
            projectedCombinedCostsForThisCombo = zeros(destinationRFinputCombinationsNum, nearbyDestinationRFinputCombinationsNum);

            % Compute projected costs
            for iDestinationRFcomboIndex = 1:destinationRFinputCombinationsNum
                % Compute the projected cost for the source RGC after the swap

                % Remove the cone(s) to be swapped out
                destinationRFinputIndicesToBeSwappedOut = theDestinationRFinputIndices(destinationRFinputCombinations(iDestinationRFcomboIndex,:));
                [destinationRFinputIndicesAfterSwap,ia] = setdiff(theDestinationRFinputIndices, destinationRFinputIndicesToBeSwappedOut); 
                destinationRFinputWeightsAfterSwap = theDestinationRFinputWeights(ia);

                projectedCombinedCostsForThisCombo(iDestinationRFcomboIndex,:) = computeCostsForAllSwapsWithNearbyDestinatonRFs(...
                    obj, ...
                    theDestinationRFindex,theNearbyDestinationRFindex, ...
                    destinationRFinputIndicesToBeSwappedOut, ...
                    destinationRFinputIndicesAfterSwap, ...
                    destinationRFinputWeightsAfterSwap, ...
                    theNearbyDestinationRFinputIndices, ...
                    theNearbyDestinationRFinputWeights, ...
                    nearbyDestinationRFinputCombinations);
            end % iDestinationRFcomboIndex
   
            % Find the (iDestinationRFcomboIndex, iNeighborComboIndex) pair that minimizes the projected costs
            [minCost,idx] = min(projectedCombinedCostsForThisCombo(:));
            [iDestinationRFcomboIndexOptimal, iNearbyComboIndexOptimal] = ind2sub(size(projectedCombinedCostsForThisCombo),idx);

            % Store the optimal input combinations for this 
            % (nearby destination RF x number of swapped inputs) as well as the corresponding cost
            optimalDestinationRFinputCombination{iNearbyDestinationRF, inputsNumSwapped} = destinationRFinputCombinations(iDestinationRFcomboIndexOptimal,:);
            optimalNearbyDestinationRFinputCombination{iNearbyDestinationRF, inputsNumSwapped} = nearbyDestinationRFinputCombinations(iNearbyComboIndexOptimal,:);
            
            % Store the optimal cost
            optimalProjectedCosts(iNearbyDestinationRF,inputsNumSwapped) = minCost;
        end % inputsNumSwapped
    end % iNearbyDestinationRF

    
    % Find the optimal nearby RGC and optimal # of swapped cones
    [costAfterOptimalSwap, idx] = min(optimalProjectedCosts(:));
    % Retrieve the optimal cone combinations for source and neighbor RGC
    [optimalNearbyDestinationRFindex, optimalInputsNumSwapped] = ind2sub(size(optimalProjectedCosts), idx);
        
    % Check to see if the cost after the optimal swap is lower than the
    % cost before swapping
    beneficialSwapWasFound = false;

   
    if (costAfterOptimalSwap < costsBeforeSwapping(optimalNearbyDestinationRFindex))
        beneficialSwapWasFound = true;

        % Retrieve the optimal input combinations for the destinationRF and the nearby destination RF
        optimalDestinationRFinputCombo = optimalDestinationRFinputCombination{optimalNearbyDestinationRFindex,optimalInputsNumSwapped};
        optimalNearbyDestinationRFinputCombo = optimalNearbyDestinationRFinputCombination{optimalNearbyDestinationRFindex,optimalInputsNumSwapped};
    
        % Retrieve the corresponding input indices
        destinationRFinputIndicesToBeSwapped = theDestinationRFinputIndices(optimalDestinationRFinputCombo);
        theNearbyDestinationRFinputIndices = allNearbyDestinationRFinputIndices{optimalNearbyDestinationRFindex};
        nearbyDestinationRFinputIndicesToBeSwapped = theNearbyDestinationRFinputIndices(optimalNearbyDestinationRFinputCombo);
        
        % Do the swap
        obj.swapInputsFromDestinationRFWithInputsOfNearbyDestinationRF(...
            destinationRFinputIndicesToBeSwapped, theDestinationRFindex, ...
            nearbyDestinationRFinputIndicesToBeSwapped, allNearbyDestinationRFindices(optimalNearbyDestinationRFindex));     
    end
end

        

function projectedCombinedCostsForThisCombo = computeCostsForAllSwapsWithNearbyDestinatonRFs(obj,...
                theDestinationRFindex,theNearbyDestinationRFindex, ...
                destinationRFinputIndicesToBeSwappedOut, ...
                destinationRFinputIndicesAfterSwap, ...
                destinationRFinputWeightsAfterSwap, ...
                theNearbyDestinationRFinputIndices, ...
                theNearbyDestinationRFinputWeights, ...
                nearbyDestinationRFinputCombinations)
       
    destinationCombinationsNum = size(nearbyDestinationRFinputCombinations,1);
    projectedCombinedCostsForThisCombo = zeros(1, destinationCombinationsNum);

    for iNeighborComboIndex = 1:destinationCombinationsNum
        % Determine inputs to be swapped out
        nearbyDestinationRFinputIndicesBeforeSwap = theNearbyDestinationRFinputIndices;
        nearbyDestinationRFinputIndicesToBeSwappedOut = nearbyDestinationRFinputIndicesBeforeSwap(nearbyDestinationRFinputCombinations(iNeighborComboIndex,:));
                    
        % Remove the inputs to be swapped out
        [nearbyDestinationRFinputIndicesAfterSwap,ia] = setdiff(nearbyDestinationRFinputIndicesBeforeSwap, nearbyDestinationRFinputIndicesToBeSwappedOut); 
        nearbyDestinationRFinputWeightsAfterSwap = theNearbyDestinationRFinputWeights;
        nearbyDestinationRFinputWeightsAfterSwap = nearbyDestinationRFinputWeightsAfterSwap(ia);

        % For the nearby destination RF: Add the cone(s) to be swapped in
        % from the destination RF
        k = numel(destinationRFinputIndicesToBeSwappedOut);
        newIndices = numel(nearbyDestinationRFinputIndicesAfterSwap)+(1:k);
        nearbyDestinationRFinputIndicesAfterSwap(newIndices) = destinationRFinputIndicesToBeSwappedOut;
        nearbyDestinationRFinputWeightsAfterSwap(newIndices) = 1;

        % For the destination RF: Add the inputs to be swapped in from the
        % nearby destination RF
        newIndices = numel(destinationRFinputIndicesAfterSwap)+(1:k);
        destinationRFinputIndicesAfterSwap(newIndices) = nearbyDestinationRFinputIndicesToBeSwappedOut;
        destinationRFinputWeightsAfterSwap(newIndices) = 1;

        % Compute projected cost for sourceRGC
        projectedCostForDestinationRF = obj.inputMaintenanceCost(...
                        destinationRFinputIndicesAfterSwap, destinationRFinputWeightsAfterSwap, ...
                        obj.destinationRFspacingsFromCentroids(theDestinationRFindex));
                    
                    
        % Compute projected cost for neighborRGC
        projectedCostForNearbyDestinationRF = obj.inputMaintenanceCost(...
                      nearbyDestinationRFinputIndicesAfterSwap, ...
                      nearbyDestinationRFinputWeightsAfterSwap, ...
                      obj.destinationRFspacingsFromCentroids(theNearbyDestinationRFindex));

        % Compute projected cost due to overlap
        projectedCostDueToRFoverlap = obj.overlappingDestinationRFCost(...
                      theNearbyDestinationRFindex, ...
                      nearbyDestinationRFinputIndicesAfterSwap, ...
                      nearbyDestinationRFinputWeightsAfterSwap, ...
                      theDestinationRFindex, ...
                      destinationRFinputIndicesAfterSwap, ...
                      destinationRFinputWeightsAfterSwap);

         % Total cost
         projectedCombinedCostsForThisCombo(1, iNeighborComboIndex) = ...
                        projectedCostForDestinationRF(1) + ...
                        projectedCostForNearbyDestinationRF(1) + ...
                        projectedCostDueToRFoverlap;

    end % iNeighborComboIndex
end

