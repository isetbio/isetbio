function optimizeTransferOfInputRFsToZeroInputDestinationRF(obj,...
             theMultiInputDestinationRFindex, theMultiInputDestinationRFinputIndices, ...
             theMultiInputDestinationRFinputWeights, theZeroInputDestinationRF)

    % Optimize how many and which of theMultiInputDestinationRFinputIndices 
    % will be transfered to theZeroInputDestinationRF
    
    % All combinations of half of the input RFs
    halfInputRFsNum = round(numel(theMultiInputDestinationRFinputIndices)/2);
    
    C = nchoosek(1:numel(theMultiInputDestinationRFinputIndices),halfInputRFsNum);
    combinationsNum = size(C,1);
    projectedCostsForAllCombinations = zeros(combinationsNum,1);

    parfor iCombo = 1:combinationsNum
        % The indices of the input RFs to be transfered
        theInputRFindices = theMultiInputDestinationRFinputIndices(C(iCombo,:));

        % Remove theInputRFindices from theMultiInputDestinationRFinputIndices
        [theRemainingInputRFindices, ia] = setdiff(theMultiInputDestinationRFinputIndices, theInputRFindices);
        theRemainingInputRFweights = theMultiInputDestinationRFinputWeights(ia);

        % Compute projected cost for theMultiInputDestinationRFindex
        theCostComponents = obj.inputMaintenanceCost(...
            theRemainingInputRFindices, theRemainingInputRFweights, ...
            obj.destinationRFspacingsFromCentroids(theMultiInputDestinationRFindex));

        projectedCostsForAllCombinations(iCombo) = theCostComponents(1);
    end

    % Find the combination that minimizes the total cost
    [~,iBestConeCombo] = min(projectedCostsForAllCombinations(:));
    
   
    % The indices of the inputRFs to be transfered
    theInputRFindices = theMultiInputDestinationRFinputIndices(C(iBestConeCombo,:));

    % Reassign the selected inputRFs from theMultiInputDestinationRFindex to theZeroInputDestinationRF

    % DISCONNECT theInputRFindices from theMultiInputDestinationRFindex
    obj.connectivityMatrix(theInputRFindices, theMultiInputDestinationRFindex) = 0;

    % And CONNECT the to theZeroInputDestinationRF
    obj.connectivityMatrix(theInputRFindices, theZeroInputDestinationRF) = 1;
    
    % Update the centroids of the 2 destination RFs
    destinationRFList = [theMultiInputDestinationRFindex theZeroInputDestinationRF];
    obj.updateDestinationCentroidsFromInputs(destinationRFList);

end

