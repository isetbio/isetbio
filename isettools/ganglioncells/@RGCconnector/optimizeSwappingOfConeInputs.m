function beneficialSwapWasFound = optimizeSwappingOfConeInputs(obj, ...
     sourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
     neighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights)

    % Compute cost for the source RGC to maintain its current inputs
    sourceRGCcost = obj.costToMaintainInputs(...
                theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
                obj.localRGCRFspacingsMicrons(sourceRGCindex));

    % Only swap up to a max of the input cones -1
    sourceRGCconeInputsNum = numel(theSourceRGCinputConeIndices);
    maxSourceRGCconesNumSwapped = max([1 sourceRGCconeInputsNum-1]);

    % Cost for all combinations of nearby RGCs and # of cones to be swapped
    optimalProjectedCosts = inf(numel(neighboringRGCindices), 100);
    optimalSourceConeCombination = cell(numel(neighboringRGCindices), 100);
    optimalNeighboringConeCombination = cell(numel(neighboringRGCindices), 100);

    % Costs before swapping
    costsBeforeSwapping = zeros(1, numel(neighboringRGCindices));
    
    for iNearbyRGC = 1:numel(neighboringRGCindices)

       % fprintf('Evaluating swaps with nearby RGC %d / %d\n',iNearbyRGC, numel(neighboringRGCindices));

        % Compute combined source + neighbor RGC cost before swapping
        neigboringRGCindex = neighboringRGCindices(iNearbyRGC);

        neighboringRGCinputConeIndices = allNeighboringRGCsInputConeIndices{iNearbyRGC};
        neighboringRGCinputConeWeights = allNeighboringRGCsInputConeWeights{iNearbyRGC};

        costsBeforeSwapping(iNearbyRGC) = sourceRGCcost + ...
                obj.costToMaintainInputs(neighboringRGCinputConeIndices, ...
                                         neighboringRGCinputConeWeights, ...
                                         obj.localRGCRFspacingsMicrons(neigboringRGCindex)) + ...
                obj.costToMaintainOverlappingInputs(...
                                         neigboringRGCindex, ...
                                         neighboringRGCinputConeIndices, ...
                                         neighboringRGCinputConeWeights, ...
                                         sourceRGCindex, ...
                                         theSourceRGCinputConeIndices, ...
                                         theSourceRGCinputConeWeights);

        % Find max # of cones to be swapped
        neigboringRGCconeInputsNum = numel(neighboringRGCinputConeIndices);
        maxNeigboringRGCconesNumSwapped = max([1 neigboringRGCconeInputsNum]);

        maxConesNumSwapped = min([...
            maxSourceRGCconesNumSwapped ...
            maxNeigboringRGCconesNumSwapped ...
            obj.wiringParams.maxNumberOfConesToSwap]);


        parfor conesNumSwapped = 1:maxConesNumSwapped
            
            % Compute all conesNumSwapped permutations for the source RGC cone inputs
            sourceConeCombinations = nchoosek(1:sourceRGCconeInputsNum, conesNumSwapped);

            % Compute all conesNumSwapped permutations for the neighboring RGC cone inputs
            neighboringConeCombinations = nchoosek(1:neigboringRGCconeInputsNum, conesNumSwapped);

            % Allocate memory to store all the projected costs
            sourceCombinationsNum = size(sourceConeCombinations,1);
            destinationCombinationsNum = size(neighboringConeCombinations,1);
            projectedCombinedCostsForThisCombo = zeros(sourceCombinationsNum, destinationCombinationsNum);

            % Compute projected costs
            for iSourceComboIndex = 1:sourceCombinationsNum
                % Compute the projected cost for the source RGC after the swap

                % Remove the cone(s) to be swapped out
                sourceRGCconeIndicesToBeSwappedOut = theSourceRGCinputConeIndices(sourceConeCombinations(iSourceComboIndex,:));
                [sourceRGCconeIndicesAfterSwap,ia] = setdiff(theSourceRGCinputConeIndices, sourceRGCconeIndicesToBeSwappedOut); 
                sourceRGCconeWeightsAfterSwap = theSourceRGCinputConeWeights(ia);

                projectedCombinedCostsForThisCombo(iSourceComboIndex,:) = evaluateAllNeighboringCombinationCosts(...
                    obj, ...
                    sourceRGCindex,neigboringRGCindex, ...
                    sourceRGCconeIndicesToBeSwappedOut, ...
                    sourceRGCconeIndicesAfterSwap, ...
                    sourceRGCconeWeightsAfterSwap, ...
                    neighboringRGCinputConeIndices, ...
                    neighboringRGCinputConeWeights, ...
                    neighboringConeCombinations);
            end % iSourceComboIndex
   
            % Find the (iSourceComboIndex, iNeighborComboIndex) pair that minimizes the projected costs
            [minCost,idx] = min(projectedCombinedCostsForThisCombo(:));
            [iSourceComboIndexOptimal, iNeighborComboIndexOptimal] = ind2sub(size(projectedCombinedCostsForThisCombo),idx);

            % Store the optimal cone combinations for this nearbyRGC and #
            % of swapped cones and the corresponding cost
            optimalSourceConeCombination{iNearbyRGC, conesNumSwapped} = sourceConeCombinations(iSourceComboIndexOptimal,:);
            optimalNeighboringConeCombination{iNearbyRGC, conesNumSwapped} = neighboringConeCombinations(iNeighborComboIndexOptimal,:);
            
            % Store the optimal cost
            optimalProjectedCosts(iNearbyRGC,conesNumSwapped) = minCost;
        end % conesNumSwapped
    end % iNearbyRGC

    
    % Find the optimal nearby RGC and optimal # of swapped cones
    [costAfterOptimalSwap, idx] = min(optimalProjectedCosts(:));
    % Retrieve the optimal cone combinations for source and neighbor RGC
    [optimalNearbyRGCindex, optimalConesNumSwapped] = ind2sub(size(optimalProjectedCosts), idx);
        
    % Check to see if the cost after the optimal swap is lower than the
    % cost before swapping
    beneficialSwapWasFound = false;

   
    if (costAfterOptimalSwap < costsBeforeSwapping(optimalNearbyRGCindex))
        beneficialSwapWasFound = true;

        % Retrieve the optimal cone combinations for source and neighbor RGC
        optimalSourceConeCombo = optimalSourceConeCombination{optimalNearbyRGCindex,optimalConesNumSwapped};
        optimalNeighboringConeCombo = optimalNeighboringConeCombination{optimalNearbyRGCindex,optimalConesNumSwapped};
    
        % Retrieve the corresponding cone indices
        sourceRGCconeIndicesToBeSwapped = theSourceRGCinputConeIndices(optimalSourceConeCombo);
        theNeighboringRGCinputConeIndices = allNeighboringRGCsInputConeIndices{optimalNearbyRGCindex};
        neighboringRGCconeIndicesToBeSwapped = theNeighboringRGCinputConeIndices(optimalNeighboringConeCombo);
        
        % Do the swap
        obj.swapConesFromSourceRGCWithConesOfDestinationRGC(...
            sourceRGCconeIndicesToBeSwapped, sourceRGCindex, ...
            neighboringRGCconeIndicesToBeSwapped, neighboringRGCindices(optimalNearbyRGCindex));     
    end


end

function projectedCombinedCostsForThisCombo = evaluateAllNeighboringCombinationCosts(obj,...
                sourceRGCindex,neigboringRGCindex, ...
                sourceRGCconeIndicesToBeSwappedOut, ...
                sourceRGCconeIndices, ...
                sourceRGCconeWeights, ...
                neighboringRGCinputConeIndices, ...
                neighboringRGCinputConeWeights, ...
                neighboringConeCombinations)
       
    destinationCombinationsNum = size(neighboringConeCombinations,1);
    projectedCombinedCostsForThisCombo = zeros(1, destinationCombinationsNum);

    for iNeighborComboIndex = 1:destinationCombinationsNum
        % Determine cones to be swapped out
        neighboringRGCconeIndicesBeforeSwap = neighboringRGCinputConeIndices;
        neighboringRGCconeIndicesToBeSwappedOut = neighboringRGCconeIndicesBeforeSwap(neighboringConeCombinations(iNeighborComboIndex,:));
                    
        % Remove the cone(s) to be swapped out
        [neighboringRGCconeIndicesAfterSwap,ia] = setdiff(neighboringRGCconeIndicesBeforeSwap, neighboringRGCconeIndicesToBeSwappedOut); 
        neighboringRGCconeWeightsAfterSwap = neighboringRGCinputConeWeights;
        neighboringRGCconeWeightsAfterSwap = neighboringRGCconeWeightsAfterSwap(ia);

        % For the neirbhbor RGC: Add the cone(s) to be swapped in from the source RGC
        k = numel(sourceRGCconeIndicesToBeSwappedOut);
        newIndices = numel(neighboringRGCconeIndicesAfterSwap)+(1:k);
        neighboringRGCconeIndicesAfterSwap(newIndices) = sourceRGCconeIndicesToBeSwappedOut;
        neighboringRGCconeWeightsAfterSwap(newIndices) = 1;

        % For the source RGC: Add the cone(s) to be swapped in from the neighbor
        sourceRGCconeIndicesAfterSwap = sourceRGCconeIndices;
        sourceRGCconeWeightsAfterSwap = sourceRGCconeWeights;
        newIndices = numel(sourceRGCconeIndicesAfterSwap)+(1:k);
        sourceRGCconeIndicesAfterSwap(newIndices) = neighboringRGCconeIndicesToBeSwappedOut;
        sourceRGCconeWeightsAfterSwap(newIndices) = 1;

        % Compute projected cost for sourceRGC
        projectedCostForSourceRGC = obj.costToMaintainInputs(...
                        sourceRGCconeIndicesAfterSwap, sourceRGCconeWeightsAfterSwap, ...
                        obj.localRGCRFspacingsMicrons(sourceRGCindex));
                    
                    
        % Compute projected cost for neighborRGC
        projectedCostForNeighboringRGC = obj.costToMaintainInputs(...
                      neighboringRGCconeIndicesAfterSwap, ...
                      neighboringRGCconeWeightsAfterSwap, ...
                      obj.localRGCRFspacingsMicrons(neigboringRGCindex));

        % Compute projected cost due to overlap
        projectedCostDueToRFoverlap = obj.costToMaintainOverlappingInputs(...
                              neigboringRGCindex, ...
                              neighboringRGCconeIndicesAfterSwap, ...
                              neighboringRGCconeWeightsAfterSwap, ...
                              sourceRGCindex, ...
                              sourceRGCconeIndicesAfterSwap, ...
                              sourceRGCconeWeightsAfterSwap);

         % Total cost
         projectedCombinedCostsForThisCombo(1, iNeighborComboIndex) = ...
                        projectedCostForSourceRGC + ...
                        projectedCostForNeighboringRGC + ...
                        projectedCostDueToRFoverlap;

    end % iNeighborComboIndex
end