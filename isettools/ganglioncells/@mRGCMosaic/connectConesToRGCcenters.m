% Method to connect cones within the FOV of the RGC mosaic to RGC RF centers
% Called by obj.wireRFcenterToInputCones()
function connectConesToRGCcenters(obj, idxConesInsideRGCmosaic, visualizeWiringStages)

    % How far away to look for nearby RGCs for donating cones
    % This is times the mean cone spacing at that location
    searchRadiusFactor = 1.2;
    
    % How far away to look for nearby orphan RGCs for donating cones
    % This is time the cell's eccentricity
    oprhanSearchRadiusFactor = 3.0;
    
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass. Also
    % some RGCs will receive zero cone inputs because the closest cone was
    % an S-cone. This method sets the 
    % - obj.coneConnectivityMatrix and returns 
    % - the distance of each RGC to its closest cone.
    connectEachConeToItsNearestRGC(obj,idxConesInsideRGCmosaic);
    
    domainLimits = [1.5 3.0 -0.4 0.4]*obj.inputConeMosaic.micronsPerDegree;
    
    if (visualizeWiringStages)
        figNo = 1;
        visualizeProgress(obj, figNo, domainLimits);
    end
    
    
    % Optimize central retina (0 to 1 degs)
    rangeDegs = [0 1];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [-1 1 -1 1]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
    
    % Optimize parafoveal retina (1 to 2 degs)
    rangeDegs = [1 2];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [1 2 -1 1]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
    
    
    % Optimize parafoveal retina (2 to 3 degs)
    rangeDegs = [2 3];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [2 3 -1 1]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
    
    
    % Optimize parafoveal retina (3 to 4 degs)
    rangeDegs = [3 4];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [3 4 -1 1]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
    
    % Optimize parafoveal retina (4 tp 5 degs)
    rangeDegs = [4 5];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [4 5 -1 1]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
   
    % Optimize peripheral retina (5 - 7 degs)
    rangeDegs = [5 7];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    visualizationDomainLimits = [5 7 -2 2]*obj.inputConeMosaic.micronsPerDegree;
    optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
        visualizeWiringStages, visualizationDomainLimits);
    
    
%     % Optimize peripheral retina
%     rangeDegs = [7 10];
%     useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
%     visualizationDomainLimits = [7 10 -3 3]*obj.inputConeMosaic.micronsPerDegree;
%     optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
%         useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
%        visualizeWiringStages, visualizationDomainLimits);
% %    
%     % Optimize peripheral retina
%     rangeDegs = [10 15];
%     useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
%     visualizationDomainLimits = [10 15 -5 5]*obj.inputConeMosaic.micronsPerDegree;
%     optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs, ...
%         useOrphanRGCsToSplitMatchedConeMultiInputRGCs, ...
%         visualizeWiringStages, visualizationDomainLimits);
    
    
end


function optimizeSubregion(obj, searchRadiusFactor, oprhanSearchRadiusFactor, rangeDegs,  ...
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeWiringStages, ...
    visualizationDomainLimits)

    fprintf(2, 'Optimizing region between %2.1f and %2.1f degrees\n', rangeDegs(1), rangeDegs(2));

    % ======================= DEALING WITH N-INPUT RGCS =======================
    
    % Pass 1. Minimize number of RGCs that connect to N cones by
    % trying to donate one cone to a nearby RGC with up to N-2  cone inputs
    
    for inputConesNum = 10:-1:4
        % Find all RGCs with inputConesNum number of input cones
        [nConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, inputConesNum);
    
        % Minimize their frequency by donating one cone to a nearby
        % RGC with up to inputsConesNum-2 input cones
        successfullReallocations = minimizeFrequencyOfNInputRGCSByDonatingOneConeToNearbyRGC(obj, ...
            searchRadiusFactor, nConeInputRGCs, coneInputIDs);
        fprintf('\t %d/%d (%d cone RGC -> %d cone RGC) reassignments.\n', ...
            successfullReallocations(1), successfullReallocations(2), ...
            inputConesNum, inputConesNum-2);
    end
    
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    % ======================= DEALING WITH 3-INPUT RGCS =======================
    
    % Pass 2A. Minimize number of RGCs that connect to 3 cones by
    % trying to donate one cone to a nearby RGC with one cone input only.
    % The nearby cone type of the single input RGC must match the donated
    % cone type
    
    % Find all RGCs with 3-cone inputs
    [threeConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, 3); 
    
    % Minimize their frequency by donating one cone to a nearby single-cone
    % RGC (with cone input matching the donated cone)
    donatedConeMustMatch = true;
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC with matched cone type) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
   
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    % Pass 2b Minimize the frequency of remaining 3-cone input RGCs by donating one cone to a nearby single-cone
    % RGC (now, without matching cone input to the donated cone)
    donatedConeMustMatch = ~true;
    % Find remaining 3-cone RGCs
    [threeConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, 3); 
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC with unmatched cone type) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
   
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    % Pass 2C. Mimize the frequency of still remaining 3-cone input RGCs by donating
    % one cone to a nearby OrphanRGC
    % Find remaining 3- cone RGCs
    [threeConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, 3); 
    successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
            threeConeInputRGCs, coneInputIDs, oprhanSearchRadiusFactor);
    fprintf('\t %d/%d (3 cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    if (visualizeWiringStages)
        figNo = 1000+100*rangeDegs(1)+1;
        visualizeProgress(obj, figNo, visualizationDomainLimits);
    end
    
    
    % ======================= DEALING WITH 2-INPUT RGCS =======================
    
    % Pass 3A.  Minimize number of RGCs that connect to two mismatched cone
    % types by re-assigning one of the cones to a nearby RGC with 1 cone input RGCs
    % of matched cone type
    
    % Find all RGCs with 2-mismatched cone inputs
    [twoMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'mismatched', rangeDegs);
    
    % Minimize their frequency
    successfullReallocations = minimizeFrequencyOf2MismatchedInputRGCSByDonatingToNearbyRGC(obj, ...
        searchRadiusFactor, twoMismatchedInputConeTypeRGCs, coneInputIDs);
    fprintf('\t %d/%d of (2 mismatched cone RGC -> nearby 1 cone RGC with matched cone) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    % Pass 3B. Mimize the frequency of still remaining 2-mismatched cone input RGCs by donating
    % one cone to a nearby OrphanRGC
    [twoMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'mismatched', rangeDegs);
    successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
            twoMismatchedInputConeTypeRGCs, coneInputIDs, oprhanSearchRadiusFactor);
    fprintf('\t %d/%d (2 mismatched cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
     
    
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    % Pass 3C. Mimize the frequency of 2-matched cone input RGCs by donating
    % one cone to a nearby OrphanRGC
    if (useOrphanRGCsToSplitMatchedConeMultiInputRGCs)
        [twoMatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'matched', rangeDegs);
        successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
            twoMatchedInputConeTypeRGCs, coneInputIDs, oprhanSearchRadiusFactor);
        fprintf('\t %d/%d (2 matched cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    end
   
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    if (visualizeWiringStages)
        figNo = 1000+100*rangeDegs(1)+2;
        visualizeProgress(obj, figNo, visualizationDomainLimits);
    end
    
    
    %
    % NEW RULE
    %
    
    % Pass 4. Look for RGCs with 2 matched cone inputs whose separation
    % is > mean local spacing, and see if we can reassign one of those
    % inputs to a nearby RGC with same cone type so that the separation of
    % these 2 cones is less 
    
    thresholdDistanceFactor = 1.1;
    [twoMatchedInputConeTypeRGCsWithLargeInputSeparation, coneInputIDs] = RGCsWithTwoConeInputsSeparatedMoreThanDesiredDistance(obj, 'matched', thresholdDistanceFactor, rangeDegs);
    
    % Minimize their frequency
    successfullReallocations = minimizeFreqOf2MatchedInputRGCSResultingInSmallerConeSeparation(obj, ...
        searchRadiusFactor, twoMatchedInputConeTypeRGCsWithLargeInputSeparation, coneInputIDs);
    fprintf('\t %d/%d of (2 matched cone RGC -> nearby 1 cone RGC with matched cone resulting in smaller cone separation) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
      
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    
    % Repeat pass 2
    % Find all RGCs with 3-cone inputs
    [threeConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, 3); 
    
    % Minimize their frequency by donating one cone to a nearby single-cone
    % RGC (with cone input matching the donated cone)
    donatedConeMustMatch = true;
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC with matched cone type) reassignments (REPEAT).\n', successfullReallocations(1), successfullReallocations(2));
   
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
    % Pass 2b Minimize the frequency of remaining 3-cone input RGCs by donating one cone to a nearby single-cone
    % RGC (now, without matching cone input to the donated cone)
    donatedConeMustMatch = ~true;
    % Find remaining 3-cone RGCs
    [threeConeInputRGCs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, 3); 
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC with unmatched cone type) reassignments (REPEAT).\n', successfullReallocations(1), successfullReallocations(2));
   
    % Update RGC positions based on their current inputs
    updateRGCpositionsBasedOnCurrentInputs(obj);
    
    
end

function updateRGCpositionsBasedOnCurrentInputs(obj)
    
    for sourceRGCIndex = 1:size(obj.coneConnectivityMatrix,2)
        % Update the position of the source RGC to be the centroid of its new cone inputs
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, sourceRGCIndex)) == 1);
        obj.rgcRFpositionsMicrons(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
        obj.rgcRFpositionsDegs(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
    end
end


function [rgcIDs, coneInputIDs] = RGCsWithTwoConeInputsSeparatedMoreThanDesiredDistance(obj, coneInputSchema, thresholdDistanceFactor, rangeDegs)   
    % List of 2-cone input RGCs
    rgcIDs = [];
    % [N x 2] indices of the 2 cone input indices to each of the 2-cone RGC
    coneInputIDs = [];

    % Find all 2-input RGCs
    twoInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 2);
    
    % Select only those that lie between [minRadiusDegs, maxRadiusDegs]
    r2 = sum((obj.rgcRFpositionsDegs(twoInputRGCindices,:)).^2,2);
    idx = find((r2>=rangeDegs(1)^2) & (r2 <= rangeDegs(2)^2));
    twoInputRGCindices = twoInputRGCindices(idx);
    
    % Determine those that have a separation > threshold
    twoInputRGCindicesSeparatedMoreThanThresholdDistance = [];
    for k = 1:numel(twoInputRGCindices)
        % Get the indices of cones connected to this RGC
        rgcIndex = twoInputRGCindices(k);
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        % Compute separation of these cones
        meanSpacing = mean(obj.inputConeMosaic.coneRFspacingsMicrons(indicesOfConeInputs));
        thresholdDistance = meanSpacing * thresholdDistanceFactor;
        actualDistance = sqrt(sum((obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs(1),:) - ...
                         obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs(2),:)).^2,2));
        if (actualDistance > thresholdDistance)
            twoInputRGCindicesSeparatedMoreThanThresholdDistance = cat(2, twoInputRGCindicesSeparatedMoreThanThresholdDistance, rgcIndex);
        end
    end
    
    for k = 1:numel(twoInputRGCindicesSeparatedMoreThanThresholdDistance)
        % Get the indices of cones connected to this RGC
        rgcIndex = twoInputRGCindicesSeparatedMoreThanThresholdDistance(k);
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        switch coneInputSchema
            case 'mismatched'
                % If cone input types differ, add to list 
                if (obj.inputConeMosaic.coneTypes(indicesOfConeInputs(1)) ~= obj.inputConeMosaic.coneTypes(indicesOfConeInputs(2)) )
                    rgcIDs = cat(2, rgcIDs, rgcIndex);
                    coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
                end
            case 'matched'
                % If cone input types agree, add to list 
                if (obj.inputConeMosaic.coneTypes(indicesOfConeInputs(1)) == obj.inputConeMosaic.coneTypes(indicesOfConeInputs(2)) )
                    rgcIDs = cat(2, rgcIDs, rgcIndex);
                    coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
                end
            otherwise
                error('Unknown coneInputSchema: ''%s''. Must be either ''matched'', or ''mismatched''.', coneInputSchema);
        end
    end
    if (isempty(coneInputIDs))
        return;
    end
    
    % Return 2-cone input RGC indices sorted according to the RGC eccentricity
    ecc = sum(obj.rgcRFpositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~=2)
        error('Size must be 2');
    end
end

function successfullReallocations = minimizeFreqOf2MatchedInputRGCSResultingInSmallerConeSeparation(obj, ...
        searchRadiusFactor, twoMatchedInputConeTypeRGCsWithLargeInputSeparation, inputConeIndices)
    
    successfullReallocations = 0;
    
    if (isempty(twoMatchedInputConeTypeRGCsWithLargeInputSeparation))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(twoMatchedInputConeTypeRGCsWithLargeInputSeparation)
        % Get the source RGC index
        sourceRGCIndex = twoMatchedInputConeTypeRGCsWithLargeInputSeparation(iRGC);
        % Ge the indices of cones connecting to it
        sourceRGCinputConeIndices = inputConeIndices(iRGC,:);
        
        
        % Select which cone to reassign.
        % Randomly select which one to try to re-assign first
        if (rand < 0.5)
            indexOfConeToBeReassigned = sourceRGCinputConeIndices(1);
            indexOfOtherConeToBeReassigned = sourceRGCinputConeIndices(2);    
        else
            indexOfConeToBeReassigned = sourceRGCinputConeIndices(2);
            indexOfOtherConeToBeReassigned = sourceRGCinputConeIndices(1);
        end
        
        % Position of the 2 input cones
        positionOfConeToBeReassigned = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfConeToBeReassigned,:);
        positionOfOtherConeToBeReassigned = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfOtherConeToBeReassigned,:);
            
        % Compute original separation of this RGC's cone inputs
        originalConeSeparation = sqrt(sum((positionOfConeToBeReassigned-positionOfOtherConeToBeReassigned).^2,2));
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(obj.inputConeMosaic.coneRFspacingsMicrons(sourceRGCinputConeIndices));
        
        % Looking for a nearby RGC with one cone input only
        maxConeInputsOfNearbyRGC = 1;
        
        % Look for neignboring RGC this type of cone input
        typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
        
        % Look for it
        destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
            indexOfConeToBeReassigned, searchRadiusMicrons);
        
        if (~isempty(destinationRGCindex))
            
            % We have a candidate nearby RGC with a single cone input
            % Compute the distance of that cone to the cone to be
            % reassigned
            indexOfSingleConeInput = find(squeeze(obj.coneConnectivityMatrix(:, destinationRGCindex)) == 1);
            if (numel(indexOfSingleConeInput) ~= 1)
                error('We have a problem with first cone input. This RGC is supposed to have a single cone input');
            end
            candicateConePosition = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfSingleConeInput,:);
            separation = sqrt(sum((positionOfConeToBeReassigned - candicateConePosition).^2,2));
            if (separation < originalConeSeparation)
                % We got it
                % Update the connectivityMatrix, by disconnecting
                %   indexOfConeToBeReassigned  FROM  rgcIndex
                % and connecting 
                %   indexOfConeToBeReassigned  to theTargetRGCindex
                updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                        indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

                % Update successfullReallocations
                successfullReallocations = successfullReallocations + 1;
            end
        else
            
            % Look for neignboring RGC this type of cone input
            typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfOtherConeToBeReassigned);
        
            % Look for it
            destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
                maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
                indexOfOtherConeToBeReassigned, searchRadiusMicrons);
        
            if (~isempty(destinationRGCindex))
                % We have a candidate nearby RGC with a single cone input
                % Compute the distance of that cone to the cone to be
                % reassigned
                indexOfSingleConeInput = find(squeeze(obj.coneConnectivityMatrix(:, destinationRGCindex)) == 1);
                if (numel(indexOfSingleConeInput) ~= 1)
                    error('We have a problem with 2nd cone input. This RGC is supposed to have a single cone input');
                end
                candicateConePosition = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfSingleConeInput,:);
                separation = sqrt(sum((positionOfOtherConeToBeReassigned - candicateConePosition).^2,2));
                
                if (separation < originalConeSeparation)
                    % We got it
                    % Update the connectivityMatrix, by disconnecting
                    %   indexOfConeToBeReassigned  FROM  rgcIndex
                    % and connecting 
                    %   indexOfConeToBeReassigned  to theTargetRGCindex
                    updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                            indexOfOtherConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

                    % Update successfullReallocations
                    successfullReallocations = successfullReallocations + 1;
                end
            end   
        end
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(twoMatchedInputConeTypeRGCsWithLargeInputSeparation)];

end


function successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
    multiInputRGCIndices, multiInputRGConeIndices, oprhanSearchRadiusFactor)

    % Find RGCs with zero cone inputs. These are RGCs for which there was
    % no cone within the threshold distance, probably because the closest cone was an S-cone.
    orphanRGCIndices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 0);
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the source RGC index
        sourceRGCIndex = multiInputRGCIndices(iRGC);
        sourceRGCinputConeIndices = multiInputRGConeIndices(iRGC,:);
            
        % Compute its radial eccentricity
        eccDegs2 = sum(obj.rgcRFpositionsDegs(sourceRGCIndex,:).^2,2);
        
        % Search radius, equal to the cell's eccentricity with a min value of 0.5 degs
        minRadiusDegs = 0.5;
        searchRadiusDegs2 = max([minRadiusDegs^2 eccDegs2*oprhanSearchRadiusFactor]);

        % Compute distances to all orphanRGCs 
        d2ToAllOrphanRGCs = sum((bsxfun(@minus, ...
            obj.rgcRFpositionsDegs(orphanRGCIndices,:), ...
            obj.rgcRFpositionsDegs(sourceRGCIndex,:))).^2,2);
            
        [d2Min, idx] = min(d2ToAllOrphanRGCs);
        if (d2Min < searchRadiusDegs2)
            % OK, found orphanRGC within the search radius. We'll use it.
            theOrphanRGCindex = orphanRGCIndices(idx);
           
            % Remove it from the list of orphanRGCindices
            orphanRGCIndices = setdiff(orphanRGCIndices, theOrphanRGCindex);
            
            % Choose which cone input to reassign 
            if (numel(sourceRGCinputConeIndices) == 2)
                % Choose the first cone (arbitrary)
                indexOfConeToBeReassigned = sourceRGCinputConeIndices(1);
            else
                donatedConeMustMatch = false;
                indexOfConeToBeReassigned = indexOfConeWithSmallestPopulation(obj, ...
                    donatedConeMustMatch, sourceRGCinputConeIndices, sourceRGCIndex);
                if (indexOfConeToBeReassigned < 0)
                    continue;
                end
            end
            

            % Update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  sourceRGCIndex
            % and connecting 
            %   indexOfConeToBeReassigned  to theOrphanRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                indexOfConeToBeReassigned, sourceRGCIndex, theOrphanRGCindex);
            
            % Update successfullReallocations
            successfullReallocations = successfullReallocations + 1;
        end
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end


function successfullReallocations = minimizeFrequencyOf2MismatchedInputRGCSByDonatingToNearbyRGC(obj, ...
        searchRadiusFactor, multiInputRGCIndices, multiInputRGConeIndices)
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the source RGC index
        sourceRGCIndex = multiInputRGCIndices(iRGC);
        % Ge the indices of cones connecting to it
        sourceRGCinputConeIndices = multiInputRGConeIndices(iRGC,:);
        
        % Select which cone to reassign. We have one L and one M cone.
        % Randomly select which one to try to re-assign first
        if (rand < 0.5)
            indexOfConeToBeReassigned = sourceRGCinputConeIndices(1);
            indexOfOtherConeToBeReassigned = sourceRGCinputConeIndices(2);
        else
            indexOfConeToBeReassigned = sourceRGCinputConeIndices(2);
            indexOfOtherConeToBeReassigned = sourceRGCinputConeIndices(1);
        end
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(obj.inputConeMosaic.coneRFspacingsMicrons(sourceRGCinputConeIndices));
        
        % Looking for a nearby RGC with one cone input only
        maxConeInputsOfNearbyRGC = 1;
        
        % Look for neignboring RGC this type of cone input
        typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
        
        % Look for it
        destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
            indexOfConeToBeReassigned, searchRadiusMicrons);
        
        if (~isempty(destinationRGCindex))
            % Update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  rgcIndex
            % and connecting 
            %   indexOfConeToBeReassigned  to theTargetRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                    indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

            % Update successfullReallocations
            successfullReallocations = successfullReallocations + 1;
        else
            
            % Look for neignboring RGC this type of cone input
            typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfOtherConeToBeReassigned);
        
            % Look for it
            destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
                maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
                indexOfOtherConeToBeReassigned, searchRadiusMicrons);
        
            if (~isempty(destinationRGCindex))
                % Update the connectivityMatrix, by disconnecting
                %   indexOfConeToBeReassigned  FROM  rgcIndex
                % and connecting 
                %   indexOfConeToBeReassigned  to theTargetRGCindex
                updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                        indexOfOtherConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

                % Update successfullReallocations
                successfullReallocations = successfullReallocations + 1;
            end   
        end
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];

end


function successfullReallocations = minimizeFrequencyOfNInputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, multiInputRGCIndices, multiInputRGConeIndices)
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    sourceRGCConesNum = size( multiInputRGConeIndices,2);
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the source RGC index
        sourceRGCIndex = multiInputRGCIndices(iRGC);
        % Ge the indices of cones connecting to it
        sourceRGCinputConeIndices = multiInputRGConeIndices(iRGC,:);
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(obj.inputConeMosaic.coneRFspacingsMicrons(sourceRGCinputConeIndices));
        
        % Go through each of the input cones and see if there is a neighboring
        % to it RGC with up to sourceRGCConesNum-2 number of cone inputs
        maxExaminedConeInputsOfNearbyRGC = sourceRGCConesNum-2;
        
        destinationRGCindex = [];
        maxConeInputsOfNearbyRGC = 0;
        while (maxConeInputsOfNearbyRGC < maxExaminedConeInputsOfNearbyRGC) && (isempty(destinationRGCindex))
            
            % Looking for a nearby RGC with this many cone inputs only
            maxConeInputsOfNearbyRGC = maxConeInputsOfNearbyRGC + 1;
            
            iCone = 0;
            while (iCone < numel(sourceRGCinputConeIndices)) && (isempty(destinationRGCindex))

                % Index of current cone
                iCone = iCone + 1;
                indexOfConeToBeReassigned = sourceRGCinputConeIndices(iCone);

                % Looking for neignboring RGC this type of cone input
                typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);

                % Look for it
                destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
                    maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
                    indexOfConeToBeReassigned, searchRadiusMicrons);

                if (isempty(destinationRGCindex))
                    
                    % Didnt find a suitable RGC. Try without matching cone type
                    typeOfConeToBeReassigned = [];
                    destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
                        maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
                        indexOfConeToBeReassigned, searchRadiusMicrons);
                
                    if (~isempty(destinationRGCindex))
                        % Update the connectivityMatrix, by disconnecting
                        %   indexOfConeToBeReassigned  FROM  rgcIndex
                        % and connecting 
                        %   indexOfConeToBeReassigned  to theTargetRGCindex
                        updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                                indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

                        % Update successfullReallocations
                        successfullReallocations = successfullReallocations + 1;
                    end
                    
                else
                    % Update the connectivityMatrix, by disconnecting
                    %   indexOfConeToBeReassigned  FROM  rgcIndex
                    % and connecting 
                    %   indexOfConeToBeReassigned  to theTargetRGCindex
                    updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                            indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

                    % Update successfullReallocations
                    successfullReallocations = successfullReallocations + 1;
                end
            end % while iCone
        end % while iNearRGCcones
        
    end % iRGC
        
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end


    
function successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        searchRadiusFactor, multiInputRGCIndices, multiInputRGConeIndices, donatedConeMustMatch)
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the source RGC index
        sourceRGCIndex = multiInputRGCIndices(iRGC);
        % Ge the indices of cones connecting to it
        sourceRGCinputConeIndices = multiInputRGConeIndices(iRGC,:);
        
        % Select which cone to reassign. We choose the one with least
        % population, so if 2L, 1M, the M cone
        [indexOfConeToBeReassigned, typeOfConeToBeReassigned] = indexOfConeWithSmallestPopulation(obj, ...
            donatedConeMustMatch, sourceRGCinputConeIndices, sourceRGCIndex);
        
        if (indexOfConeToBeReassigned<0)
            % This means that the cone selected would be at the RF center,
            % so dont reassign it as this would create a hole in the RF
            % Eliminate that cone from the list of possible cones
            sourceRGCinputConeIndices = setdiff(sourceRGCinputConeIndices, -indexOfConeToBeReassigned);
            % And try again
            [indexOfConeToBeReassigned, typeOfConeToBeReassigned] = indexOfConeWithSmallestPopulation(obj, ...
                donatedConeMustMatch, sourceRGCinputConeIndices, sourceRGCIndex);
            if (indexOfConeToBeReassigned<0)
                % Eliminate that cone from the list of possible cones
                sourceRGCinputConeIndices = setdiff(sourceRGCinputConeIndices, -indexOfConeToBeReassigned);
                % And try with the final cone
                indexOfConeToBeReassigned = sourceRGCinputConeIndices(1);
                typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
            end
        end
        
        % If we do not care about matching the type of the donated cone, make its type []
        if (donatedConeMustMatch == false)
            typeOfConeToBeReassigned = [];
        end
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(obj.inputConeMosaic.coneRFspacingsMicrons(sourceRGCinputConeIndices));
        
        % Looking for a nearby RGC with one cone input only
        maxConeInputsOfNearbyRGC = 1;
        
        % Retrieve neighboring RGCs with no more than these cone inputs: maxConeInputsOfNearbyRGC,
        % with the majority of these inputs being of cone type: typeOfConeToBeReassigned
        destinationRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
            indexOfConeToBeReassigned, searchRadiusMicrons);
            
        if (isempty(destinationRGCindex)) 
            continue;
        end
        
        % Update the connectivityMatrix, by disconnecting
        %   indexOfConeToBeReassigned  FROM  rgcIndex
        % and connecting 
        %   indexOfConeToBeReassigned  to theTargetRGCindex
        updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);
            
        % Update successfullReallocations
        successfullReallocations = successfullReallocations + 1;
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end


function updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex)
    % DISCONNECT cone from its RGC
    if (obj.coneConnectivityMatrix(indexOfConeToBeReassigned, sourceRGCIndex) == 1)
        obj.coneConnectivityMatrix(indexOfConeToBeReassigned, sourceRGCIndex) = 0; % disconnect
    else
        error('Cone %d was not connected to RGC %d\n', indexOfConeToBeReassigned, sourceRGCIndex);
    end
    
    % And CONNECT it to the new RGC
    obj.coneConnectivityMatrix(indexOfConeToBeReassigned, destinationRGCindex) = 1;
    
    % Update the position of the source RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, sourceRGCIndex)) == 1);
    obj.rgcRFpositionsMicrons(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
            
    % Update the position of the destination RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, destinationRGCindex)) == 1);
    obj.rgcRFpositionsMicrons(destinationRGCindex,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(destinationRGCindex,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
end

% Retrieve neighboring to target cone RGCs with no more than this many cone inputs: maxConeInputsOfNearbyRGC,
% with the majority of these inputs being of cone type: typeOfConeToBeReassigned
function targetRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
       maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
       indexOfConeToBeReassigned, searchRadiusMicrons)
   
    sourceConePositionMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfConeToBeReassigned,:);
    searchRadiusForNearbyRGCs = 1.3*searchRadiusMicrons;
    d2 = sum((bsxfun(@minus, obj.rgcRFpositionsMicrons, sourceConePositionMicrons)).^2,2);
    indicesOfRGCsWithinSearchRadius = find((d2 <= searchRadiusForNearbyRGCs^2)&(d2 > 0));
    
    if (numel(indicesOfRGCsWithinSearchRadius) > 1)
        % Sort them in ascending distance order
        [~,idx] = sort(d2(indicesOfRGCsWithinSearchRadius), 'ascend');
        indicesOfRGCsWithinSearchRadius = indicesOfRGCsWithinSearchRadius(idx);
    end
    

    targetRGCindices = [];
    targetDistances = [];
    
    % Go through all of these RGCs and measure the distances of their cone
    % inputs to the cone that is to be reassigned
    for iRGC = 1:numel(indicesOfRGCsWithinSearchRadius)
        
        % Get the rgc index
        rgcIndex = indicesOfRGCsWithinSearchRadius(iRGC);
        
        % Get the indices of cones connected to this RGC
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        
        % Determine if this RGC will meet the matching requirements 
        if (numel(indicesOfConeInputs) <= maxConeInputsOfNearbyRGC) && (numel(indicesOfConeInputs)>0)
            
            % Find the closest of these cones to the source cone
            d2 = sum((bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:), sourceConePositionMicrons)).^2,2);
            [minDistance, idx] = min(d2);
            
            c1 = isempty(typeOfConeToBeReassigned);
            c2 = obj.inputConeMosaic.coneTypes(indicesOfConeInputs(idx)) == typeOfConeToBeReassigned;

            if (c1 || (~c1 && c2))
                % Add this to RGC and cone to the pool
                targetRGCindices = cat(2, targetRGCindices, rgcIndex);
                targetDistances = cat(2, targetDistances, minDistance);
            end
        end
    end
    
    % Pick the RGC with the min cone distance
    if (~isempty(targetRGCindices))
        [~,idx] = min(targetDistances);
        targetRGCindex = targetRGCindices(idx);
    else
        targetRGCindex = [];
    end
    
end



function [indexOfConeToBeReassigned, typeOfConeToBeReassigned] = indexOfConeWithSmallestPopulation(obj, ...
    donatedConeMustMatch, mismatchedConeIndices, parentRGCindex)
    % Choose the cone with the smallest population, so if we have 2 Lcones
    % and 1 M cone, choose the M cone. But only if the 3 cones are
    % equidistant from the RGC RF center. If they are not, make sure we
    % do not pick the one that is at the RF center, as that would leave an
    % RF with a hole in its center

    if (donatedConeMustMatch)
        lconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.LCONE_ID);
        mconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.MCONE_ID);
        if (numel(lconeIndices) == 0)
            indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(mconeIndices), parentRGCindex);
            assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                'CHECK 1: cone index is messed up');
            typeOfConeToBeReassigned = cMosaic.MCONE_ID;
        elseif (numel(mconeIndices) == 0)
            indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(lconeIndices), parentRGCindex);
            typeOfConeToBeReassigned = cMosaic.LCONE_ID;
            assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                'CHECK 2: cone index is messed up');
        elseif (numel(lconeIndices) < numel(mconeIndices))
            if (numel(lconeIndices) == 1)
                % See if this L-cone is closest to the RF center. If it is,
                % do not accept it, as it may create a whole
                targetConeDistanceFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(lconeIndices(1)),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                nonTargetConeDistancesFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(mconeIndices),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                if (targetConeDistanceFromRFcenter < 0.7*mean(nonTargetConeDistancesFromRFcenter))
                    % dont accept, too close to RF center
                    indexOfConeToBeReassigned = -mismatchedConeIndices(lconeIndices(1));
                else
                    % accept it
                    indexOfConeToBeReassigned = mismatchedConeIndices(lconeIndices(1));
                    assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                    'CHECK 3: cone index is messed up');
                end
            else
                indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(lconeIndices), parentRGCindex);
                assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                    'CHECK 4: cone index is messed up');
            end
            typeOfConeToBeReassigned = cMosaic.LCONE_ID;
        else
            if (numel(mconeIndices) == 1)
                % See if this M-cone is closest to the RF center. If it is,
                % do not accept it as it may create a whole
                targetConeDistanceFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(mconeIndices(1)),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                nonTargetConeDistancesFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(lconeIndices),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                if (targetConeDistanceFromRFcenter < 0.7*mean(nonTargetConeDistancesFromRFcenter))
                    % dont accept, too close to RF center
                    indexOfConeToBeReassigned = -mismatchedConeIndices(mconeIndices(1));
                else
                    % accept first one
                    indexOfConeToBeReassigned = mismatchedConeIndices(mconeIndices(1));
                    assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                       'CHECK 5: cone index is messed up');
                end
            else
                indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(mconeIndices), parentRGCindex);
                assert(obj.coneConnectivityMatrix(indexOfConeToBeReassigned, parentRGCindex) == 1, ...
                   'CHECK 6: cone index is messed up');
            end
            typeOfConeToBeReassigned = cMosaic.MCONE_ID;
        end
    else
        indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices, parentRGCindex);
        typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
    end
end


function coneIndex = mostRemoteConeFromRFcenter(obj, candidateConeIndices, parentRGCindex)
    % Find the distances of these candidate cones from the RF center
    distancesFromRFcenter = sum((bsxfun(@minus, ...
            obj.inputConeMosaic.coneRFpositionsMicrons(candidateConeIndices,:), ...
            obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
    % Find cone that is most remote to the RF center
    [~,idx] = max(distancesFromRFcenter);
    coneIndex = candidateConeIndices(idx);
end


function [rgcIDs, coneInputIDs] = RGCsWithTwoConeInputs(obj, coneInputSchema, rangeDegs)   
    % List of 2-cone input RGCs
    rgcIDs = [];
    % [N x 2] indices of the 2 cone input indices to each of the 2-cone RGC
    coneInputIDs = [];

    % Find all 2-input RGCs
    twoInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 2);
    
    % Select only those that lie between [minRadiusDegs, maxRadiusDegs]
    r2 = sum((obj.rgcRFpositionsDegs(twoInputRGCindices,:)).^2,2);
    idx = find((r2>=rangeDegs(1)^2) & (r2 <= rangeDegs(2)^2));
    twoInputRGCindices = twoInputRGCindices(idx);
    
    for k = 1:numel(twoInputRGCindices)
        % Get the indices of cones connected to this RGC
        rgcIndex = twoInputRGCindices(k);
        % Get the indices of cones connected to this RGC
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        switch coneInputSchema
            case 'mismatched'
                % If cone input types differ, add to list 
                if (obj.inputConeMosaic.coneTypes(indicesOfConeInputs(1)) ~= obj.inputConeMosaic.coneTypes(indicesOfConeInputs(2)) )
                    rgcIDs = cat(2, rgcIDs, rgcIndex);
                    coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
                end
            case 'matched'
                % If cone input types agree, add to list 
                if (obj.inputConeMosaic.coneTypes(indicesOfConeInputs(1)) == obj.inputConeMosaic.coneTypes(indicesOfConeInputs(2)) )
                    rgcIDs = cat(2, rgcIDs, rgcIndex);
                    coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
                end
            otherwise
                error('Unknown coneInputSchema: ''%s''. Must be either ''matched'', or ''mismatched''.', coneInputSchema);
        end
    end
    if (isempty(coneInputIDs))
        return;
    end
    
    % Return 2-cone input RGC indices sorted according to the RGC eccentricity
    ecc = sum(obj.rgcRFpositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~=2)
        error('Size must be 2');
    end
end


function [rgcIDs, coneInputIDs] = RGCsWithNConeInputs(obj, rangeDegs, nConeInputsNum)   
    % List of nConeInputsNum -cone input RGCs
    rgcIDs = [];
    % [N x nConeInputsNum] indices of the nConeInputsNum cone input indices to each of the nConeInputsNum RGC
    coneInputIDs = [];

    % Find all nConeInputsNum input RGCs
    nInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == nConeInputsNum);
    
    % Select only those that lie between [minRadiusDegs, maxRadiusDegs]
    r2 = sum((obj.rgcRFpositionsDegs(nInputRGCindices,:)).^2,2);
    idx = find((r2>=rangeDegs(1)^2) & (r2 <= rangeDegs(2)^2));
    nInputRGCindices = nInputRGCindices(idx);
    
    for k = 1:numel(nInputRGCindices)
        rgcIndex = nInputRGCindices(k);
        % Get the indices of cones connected to this RGC
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        % Add to list 
        rgcIDs = cat(2, rgcIDs, rgcIndex);
        coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
    end
    if (isempty(coneInputIDs))
        return;
    end
    
    % Return input RGC indices sorted according to the RGC eccentricity
    ecc = sum(obj.rgcRFpositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~= nConeInputsNum)
        error('Size must be %d', nConeInputsNum);
    end
end


function visualizeProgress(obj, figNo, domainLimits)

    hFig = figure(figNo);
    axesHandle = [];
    domain = 'microns';

    obj.visualizeConeMosaicTesselation( ...
        obj.inputConeMosaic.coneRFpositionsMicrons, obj.inputConeMosaic.coneRFspacingsMicrons, ...
        obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, domain, ...
        'figureHandle', hFig, ...
        'axesHandle', axesHandle, ...
        'showConnectedCones', true, ...
        'visualizationLimits', domainLimits, ...
        'plotTitle', 'after alignment');


    % Report connectivity stats
    obj.connectivityStats(1000 + figNo);

end


function connectEachConeToItsNearestRGC(obj, idxConesInsideRGCmosaic)

    % We only do this for cones within the FOV of the mosaic
    % (idxConesInsideRGCmosaic), not all the cones in the mosaic.
    conesNum = numel(idxConesInsideRGCmosaic);
    rgcsNum = size(obj.rgcRFpositionsMicrons,1);
    coneInputsNum = zeros(1,rgcsNum);
    
    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestRGCindices = [];
    nonSconeIndices = [];
    
    for iCone = 1:conesNum
        % Cone index in the cMosaic
        coneIndex = idxConesInsideRGCmosaic(iCone);
        
        % Do not connect S-cones to midget RGC centers
        if (obj.inputConeMosaic.coneTypes(coneIndex) == cMosaic.SCONE_ID)
            continue;
        end
        
        % Find the index of the closest RGC and connect the cone to it
        [d2, nearestRGCIndex] = min(sum((bsxfun(@minus, obj.rgcRFpositionsMicrons, obj.inputConeMosaic.coneRFpositionsMicrons(coneIndex,:)).^2),2));
       
        d = sqrt(d2);
        maxDistanceBetweenConeAndRGC = 3.0*obj.rgcRFspacingsMicrons(nearestRGCIndex);
        if (d > maxDistanceBetweenConeAndRGC)
            %fprintf('Cone is too far from nearest RGC. Will not get connected to any RGC.\n');
            continue
        end
        
        % Accumulate indices for sparse array construction 
        nonSconeIndices = cat(2, nonSconeIndices, coneIndex);
        nearestRGCindices = cat(2, nearestRGCindices, nearestRGCIndex);
        
        % Update number of cones connected to this RGC
        coneInputsNum(nearestRGCIndex) = coneInputsNum(nearestRGCIndex)+1;
    end % iCone
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    conesNumFullMosaic = size(obj.inputConeMosaic.coneRFpositionsMicrons,1);
    obj.coneConnectivityMatrix = sparse(...
        nonSconeIndices, nearestRGCindices, ones([1 numel(nonSconeIndices)]), conesNumFullMosaic, rgcsNum);
    
end

