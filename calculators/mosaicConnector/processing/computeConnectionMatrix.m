function [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
    computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, ...
    RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
    coneTypes, desiredConesToRGCratios,  orphanRGCpolicy, maximizeConeSpecificity, visualizeProcess)
        
    % Step1. Align each RGC with its nearest cone. This ensure all RGC's
    % are connected to at least one cone. Since cones are more numerous
    % than RGCs some cones will not connect to an RGC at this step. This 
    % step occurs only for RGCs for which the cone-to-RGC ratio is [1..2]
    
    RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        coneTypes, desiredConesToRGCratios, visualizeProcess);
    
    % Step 2. Connect L and M cones to midget RGC centers.
    [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = connectConesToMidgetRGCRFcenters(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        orphanRGCpolicy, maximizeConeSpecificity, ...
        coneTypes, desiredConesToRGCratios, visualizeProcess);
end

function  RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
    conePositionsMicrons, coneSpacingsMicrons, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    coneTypes, desiredConesToRGCratios, visualizeProcess)
   
    phaseString = 'Aligning: ';
    
    if (visualizeProcess)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
        
        visualizeEachAlignment = ~true;
        if (visualizeEachAlignment)
            xOutline = cosd(0:10:360);
            yOutline = sind(0:10:360);
        end
    end
    
    
    % Align RGC with cones only for those RGCs whose desired cone-to-RGC
    % ratio is less than 2
    indicesOfRGCsrequiringAlignment = find(desiredConesToRGCratios < 2); 
    rgcsNum = size(RGCRFPositionsMicrons,1);
    fprintf('%s Will align %d of %d RGCs in this patch (those RGCs have a cone-to-RGC ratio < 2)\n', ...
        phaseString, numel(indicesOfRGCsrequiringAlignment), rgcsNum);

    % Sort according to ecc
    ecc = sqrt(sum(RGCRFPositionsMicrons(indicesOfRGCsrequiringAlignment,:).^2,2));
    [~,idx] = sort(ecc, 'ascend');
    indicesOfRGCsrequiringAlignment = indicesOfRGCsrequiringAlignment(idx);
    
    % We only align mRGCs to L/M cones, not S
    eligibleConeIndices = find(coneTypes < 4);
    conePositionsMicrons = conePositionsMicrons(eligibleConeIndices,:);
    
    % Numbers of eligibe cones
    conesNum = size(conePositionsMicrons,1);
    
    % Keep a track of cones that have already been aligned to some RGC
    coneAlignedWithRGCalready = false(1, conesNum);
      
    % Go through all the RGCs one by one
    notAlignedRGCs = 0;
    for iRGC = 1:numel(indicesOfRGCsrequiringAlignment)
        % Get RGC index
        rgcIndex = indicesOfRGCsrequiringAlignment(iRGC);

        % Compute distance of this RGC to all the cones
        rgcPMicrons = RGCRFPositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPMicrons ).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 0.8*(RGCRFSpacingsMicrons(rgcIndex)));
        
        % Sort neigboring cones according to their distance to the RGC
        [~, sortedIndices] = sort(distances(coneIndicesWithinReach));
        
        % Find which cone (within the above neigborhhod) to align to. This
        % will be the closest cone that has not already been aligned with
        % another RGC
        keepGoing = true; k = 0;
        alignmentConeIndex = nan;
        while (keepGoing) && (k < numel(sortedIndices))
            k = k + 1;
            theConeIndex = coneIndicesWithinReach(sortedIndices(k));
            if (~coneAlignedWithRGCalready(theConeIndex))
                keepGoing = false;
                alignmentConeIndex = theConeIndex;
                coneAlignedWithRGCalready(theConeIndex) = true;
            end
        end
        
        % If all neighboring cones are already aligned to other RGCs, keep
        % original position
        if (isnan(alignmentConeIndex))
            notAlignedRGCs = notAlignedRGCs+1;
            %fprintf('%s Could not find a cone to align RGC %d that has not been aligned with another RGC already. Keeping original position %2.1f,%2.1f)\n', phaseString, rgcIndex, rgcPMicrons(1), rgcPMicrons(2));
            continue;
        end
        
        % Position of the cone to be aligned to
        rgcPMicronsConeAligned = conePositionsMicrons(alignmentConeIndex,:); 
        
        % Update position of RGC to get closer to  rgcPMicronsConeAligned
        % depending on desiredConesToRGCratios(rgcIndex);
        g = min([1 desiredConesToRGCratios(rgcIndex)-1]);
        rgcPMicronsAligned = g * rgcPMicrons + (1-g)*rgcPMicronsConeAligned;
        
        % Update RGC mosaic
        RGCRFPositionMicronsBeforeAlignment = RGCRFPositionsMicrons(rgcIndex,:);
        RGCRFPositionsMicrons(rgcIndex,:) = rgcPMicronsAligned;

        if (visualizeProcess)
            if (visualizeEachAlignment)
                visualizeCurrentRGCalignment(rgcPMicrons, rgcPMicronsAligned, RGCRFSpacingsMicrons(rgcIndex), ...
                    xOutline, yOutline, coneIndicesWithinReach, conePositionsMicrons, coneSpacingsMicrons, desiredConesToRGCratios);
            end
            
            % Update alignment vectors
            X1 = cat(2, X1, RGCRFPositionMicronsBeforeAlignment(1));
            X2 = cat(2, X2, rgcPMicronsAligned(1));
            Y1 = cat(2, Y1, RGCRFPositionMicronsBeforeAlignment(2));
            Y2 = cat(2, Y2, rgcPMicronsAligned(2));
            
            if (iRGC == numel(indicesOfRGCsrequiringAlignment))
                visualizeFinalAlignment(conePositionsMicrons, RGCRFPositionsMicrons, X1, X2, Y1, Y2, desiredConesToRGCratios);
            end
            
        end % visualizeProcess
    end % iRGC 
    
    fprintf('%s Could not find a cone that was not already aligned to another RGC for %d out of %d RGCs\n', phaseString, notAlignedRGCs, numel(indicesOfRGCsrequiringAlignment));
           
end

function [connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
    connectConesToMidgetRGCRFcenters(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        orphanRGCpolicy, maximizeConeSpecificity, coneTypes, desiredConesToRGCratios, visualizeProcess)
    
    
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    [connectionMatrix, numberOfConeInputs] = ...
        performPass1(conePositionsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons);
    
    % Plot ecc of orphanRGCs
    %plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, 1);
    
    
    % Second pass. For RGCs with several (> 4) cone inputs, see if we can
    % assign some of the inputs to nearby RGCs with less cone inputs
    [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass2(conePositionsMicrons,  RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs);

    conesNum = size(connectionMatrix,1);
    
    % Plot ecc of orphanRGCs
    %plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, 2);
    
    % Third pass - Maximize the frequency by which cones to 2-input RGCs have matched types
    [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass3(conePositionsMicrons,coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, maximizeConeSpecificity);
    
    % Plot ecc of orphanRGCs
    %plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, 3);
    
    % Fourth pass - In 3-input RGCs, see if we can re-assign 1 of
    % these inputs to an nearby RGC with only 1 cone input, ensuring that the cone types match
    [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass4(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, maximizeConeSpecificity);
    
    % Plot ecc of orphanRGCs
    %plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, 4);
    
    % Final phase. For RGCs with several (> 5) cone inputs, see if we can
    % assign some of the most distant inputs to nearby RGCs with less cone inputs
    [connectionMatrix, numberOfConeInputs] = ...
        performPass5(conePositionsMicrons,  RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs);
    
    % Check that each and every (non S-cone) connects to one only mRGC
    % S-cones do not connect to mRGCs
    for iCone = 1:conesNum
        checkConnectionMatrix(connectionMatrix, coneTypes, conePositionsMicrons, iCone, 'after pass 5')
    end
    
    % Plot ecc of orphanRGCs
    %plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, 5);
    
    % Remove all remaining orphan RGCs
    orphanRGCindices = find(numberOfConeInputs == 0);
    orphanRGCindicesFromConnectionMatrix = find(squeeze(sum(connectionMatrix,1)) == 0);
    if (numel(orphanRGCindicesFromConnectionMatrix) ~= numel(orphanRGCindices))
        error('Number of RFs with zero inputs do not agree');
    end
    
    for k = 1:numel(orphanRGCindicesFromConnectionMatrix)
        if (~ismember(orphanRGCindicesFromConnectionMatrix(k), orphanRGCindices))
            error('RFs with zero inputs do not agree');
        end
    end
    
    fprintf('Orphan RGCs to be removed: %d\n', numel(orphanRGCindices));
    RGCindicesToKeep = setdiff(1:size(RGCRFPositionsMicrons,1), orphanRGCindices);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(RGCindicesToKeep,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(RGCindicesToKeep);
    connectionMatrix = connectionMatrix(:,RGCindicesToKeep);
            
    fprintf('Checking connectivity matrix\n');
    parfor iCone = 1:conesNum
        checkConnectionMatrix(connectionMatrix, coneTypes, conePositionsMicrons, iCone, 'at end')
    end
    
end

function [connectionMatrix, numberOfConeInputs] = performPass1(...
    conePositionsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons)
    
    global SCONE_ID
    
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    numberOfConeInputs = zeros(1,rgcsNum);
    
    fprintf('\n-PASS 1: connecting each cone to their closest RGC ...');
    tic
    
    % Indices for sparse arrays
    closestRGCindices = [];
    nonSconeIndices = [];
    distances = [];
    
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    for iCone = 1:conesNum
        % Do not connect S-cones to RGC centers
        if (coneTypes(iCone) == SCONE_ID)
            continue;
        end
        
        % Compute distance of this cone to all RGCs
        conePosMicrons = conePositionsMicrons(iCone,:);
        % Find the index of the closest RGC and connect the iCone to it
        [d, closestRGCIndex] = min(sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2)));
       
        maxDistanceBetweenConeAndRGC = 3.0*RGCRFSpacingsMicrons(closestRGCIndex);
        if (d > maxDistanceBetweenConeAndRGC)
            fprintf('Cone is too far from nearest RGC. Will not get connected to any RGC.\n');
            continue
        end
        
        % Accumulate indices for sparse array construction 
        nonSconeIndices = cat(2, nonSconeIndices, iCone);
        closestRGCindices = cat(2, closestRGCindices, closestRGCIndex);
        distances = cat(2, distances, d);
        
        % Count number of cones in each RGC
        numberOfConeInputs(closestRGCIndex) = numberOfConeInputs(closestRGCIndex)+1;
    end % for iCone
    
    % Generate sparse connection and distance matrices
    connectionMatrix = sparse(nonSconeIndices, closestRGCindices, ones([1 numel(nonSconeIndices)]), conesNum, rgcsNum);  % connections
    
    fprintf('\nCompleted successfully in %f minutes\n',  toc/60);
end

function [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass2(conePositionsMicrons,  RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs)
    
    % Second pass. For RGCs with several (> 4) cone inputs, see if we can
    % assign some of the inputs to nearby RGCs with less cone inputs
    
    minConeInputsNum = 4;
    maxConeInputsNum = Inf;
    [rgcIDsWithSeveralConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(...
        connectionMatrix, minConeInputsNum, maxConeInputsNum, RGCRFPositionsMicrons);
    
    fprintf('\n -PASS 2: Reassigning cone inputs of %d multi-input RGCs (>4 cone inputs) to nearby RGCs with fewer cone inputs ...', numel(rgcIDsWithSeveralConeInputs));
    tic
    
    for iRGC = 1:numel(rgcIDsWithSeveralConeInputs)
        % Get the RGC index and its cone input indices
        multiInputRGCindex = rgcIDsWithSeveralConeInputs(iRGC);
        multiInputRGCconeIndices = indicesOfConnectedCones{iRGC};

        % Retrieve all the neigboring RGCs and their cone input IDs
        [neighboringRGCindices, neighboringRGCconeInputIndices] = ...
            findNeighboringRGCsAndTheirConeInputs(multiInputRGCindex, ...
            1.4*RGCRFSpacingsMicrons(multiInputRGCindex), ...
            RGCRFPositionsMicrons, connectionMatrix);
        
        meanInputsNumToNeighbors = 0;
        nonOrphanNeighbors = 0;
        for k = 1:numel(neighboringRGCindices)
            neighboringRGCconeInputIDs = neighboringRGCconeInputIndices{k};
            inputsNum = numel(neighboringRGCconeInputIDs);
            if (inputsNum>0)
                meanInputsNumToNeighbors = meanInputsNumToNeighbors + inputsNum;
                nonOrphanNeighbors = nonOrphanNeighbors + 1;
            end
        end
        meanInputsNumToNeighbors = meanInputsNumToNeighbors/nonOrphanNeighbors;
        
        for k = 1:numel(neighboringRGCindices)
            neighboringRGCindex = neighboringRGCindices(k);
            neighboringRGCconeInputIDs = neighboringRGCconeInputIndices{k};
            
            if ((numel(neighboringRGCconeInputIDs)+1 < numel(multiInputRGCconeIndices)) && ...
                (meanInputsNumToNeighbors < numel(multiInputRGCconeIndices))) 
                % Neighbor has less cone inputs and we still have more cone
                % inputs than our neighbors on average
                
                % Determine which cone input to give to this neighbor
                if (numel(neighboringRGCconeInputIDs) == 0)
                    % Assigning to an orphanRGC, so just pick the first cone
                    indexOfConeToBeReassigned = multiInputRGCconeIndices(1);
                    % update the position of the previously orphan RGC
                    RGCRFPositionsMicrons(neighboringRGCindex,:) = conePositionsMicrons(indexOfConeToBeReassigned,:);
                else
                    %  fprintf('Assigning one input for cell #%d which had %d inputs to a nearby cell with %d inputs\n', ...
                    %      multiInputRGCindex, numel(multiInputRGCconeIndices), numel(neighboringRGCconeInputIDs));
                    % Assigning to a non-orphan RGC, so pick the cone
                    % closest to the cones of that RGC
                    distances = zeros(1,numel(multiInputRGCconeIndices));
                    for kk = 1:numel(multiInputRGCconeIndices)
                        distances(kk) = min(sqrt(sum((bsxfun(@minus, conePositionsMicrons(neighboringRGCconeInputIDs,:), conePositionsMicrons(multiInputRGCconeIndices(kk),:))).^2,2)));
                    end % kk
                    [~,idx] = min(distances);
                    indexOfConeToBeReassigned = multiInputRGCconeIndices(idx);
                end
                
                % Update the cone indices to the multiInputRGC
                multiInputRGCconeIndices = setdiff(multiInputRGCconeIndices, indexOfConeToBeReassigned);
                
                % disconnect indexOfConeToBeReassigned from the multiInputRGCindex
                % and connect it to the neighboringRGCindex
                [connectionMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                    indexOfConeToBeReassigned, multiInputRGCindex, ...
                    indexOfConeToBeReassigned, neighboringRGCindex, ...
                    connectionMatrix,  numberOfConeInputs); 
            end
        end % for k
    end % iRGC
    
    fprintf('\nCompleted successfully for %d multi (>=4) input RGCs in %f minutes\n', numel(rgcIDsWithSeveralConeInputs), toc/60);
end
    
function [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass3(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, maximizeConeSpecificity)
    
    [rgcIDsWithTwoMismatchedConeInputs, indicesOfMismatchedCones] = findRGCsWithTwoMismatchedConeInputs(connectionMatrix, RGCRFPositionsMicrons,coneTypes);
    if (isempty(indicesOfMismatchedCones))
        fprintf('\n -PASS 3: skipped - no RGCs with 2 mismatched cone inputs ...');
        return;
    end
    
    fprintf('\n -PASS 3: reassigning cones in %d RGCs with 2 mismatched cone inputs to neighboring RGCs with 1 or 0 cone inputs ...',  numel(rgcIDsWithTwoMismatchedConeInputs));
    tic
    
    successfullReassignementsNum = 0;
    orphanAssignmentsNum = 0;
    
    eccDegs =  WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCRFPositionsMicrons(rgcIDsWithTwoMismatchedConeInputs,:).^2,2))/1000.0);
    
    % Probabilities of cone-specific wiring
    pAttemptConeSpecificity = rand(1,numel(rgcIDsWithTwoMismatchedConeInputs))<=maximizeConeSpecificity/100;
    
    for iRGC = 1:numel(rgcIDsWithTwoMismatchedConeInputs)
        
         % Check whether to attempt cone-specific wiring in the center
        if (pAttemptConeSpecificity(iRGC) == 0)
            continue;
        end
        
        % Get the RGC index
        mismatchedConeInputRGCindex = rgcIDsWithTwoMismatchedConeInputs(iRGC);
        mismatchedConeIndices = indicesOfMismatchedCones(iRGC,:);
        
        if (eccDegs(iRGC) <= 1.0)
            multiplier = 41;
        else
            multiplier = 7;
        end
        % Seek to connect one cone input to a nearby orphan RGC
        % Since these are orphans look in a large region (7xlocal spacing)
        orphanRGCindex = findNeigboringOrphanRGC(mismatchedConeInputRGCindex, ...
            multiplier*RGCRFSpacingsMicrons(mismatchedConeInputRGCindex), RGCRFPositionsMicrons, connectionMatrix);
       
        if (~isempty(orphanRGCindex))
            orphanAssignmentsNum = orphanAssignmentsNum + 1;
            % Reassign the first cone input (arbitary)
            indexOfConeToBeReassigned = mismatchedConeIndices(1);
            % disconnect indexOfConeToBeReassigned from the mismatchedConeInputRGCindex
            % and connect it to the orphanRGCindex
            [connectionMatrix,  numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                indexOfConeToBeReassigned, mismatchedConeInputRGCindex, ...
                indexOfConeToBeReassigned, orphanRGCindex, ...
                connectionMatrix, numberOfConeInputs);
            % update the position of the previously orphan RGC
            RGCRFPositionsMicrons(orphanRGCindex,:) = conePositionsMicrons(indexOfConeToBeReassigned,:);
        else
            % Second chance: See if we can re-assign one of the 2 mismatched cone inputs to a nearby RGC with 1 cone input
            foundRGCwithMatchedCone = false;
            for coneInputNo = 1:2
                if (foundRGCwithMatchedCone)
                    continue;
                end

                indexOfConeToBeReassigned = mismatchedConeIndices(coneInputNo);
                typeOfConeToBeReassigned = coneTypes(indexOfConeToBeReassigned);

                matchedConeInputRGCindex = findNeigboringRGCWithOneConeInputOfMatchedType(...
                    indexOfConeToBeReassigned, typeOfConeToBeReassigned, ...
                    1.5*coneSpacingsMicrons(indexOfConeToBeReassigned), ...
                    coneTypes, conePositionsMicrons, connectionMatrix);

                if (~isempty(matchedConeInputRGCindex))
                    foundRGCwithMatchedCone = true;
                    successfullReassignementsNum = successfullReassignementsNum + 1;
                    % disconnect indexOfConeToBeReassigned from the mismatchedConeInputRGCindex
                    % and connect it to the matchedRGCindex
                    [connectionMatrix,  numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                        indexOfConeToBeReassigned, mismatchedConeInputRGCindex, ...
                        indexOfConeToBeReassigned, matchedConeInputRGCindex, ...
                        connectionMatrix, numberOfConeInputs); 
                end
            end % for coneInputNo
        end
    end %iRGC
    
    fprintf('\nCompleted successfully for %d RGCs + %d orphan RGCs in %f minutes\n', successfullReassignementsNum, orphanAssignmentsNum, toc/60);

end


function [connectionMatrix, numberOfConeInputs, RGCRFPositionsMicrons] = ...
        performPass4(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, maximizeConeSpecificity)
    
    global LCONE_ID
    global MCONE_ID
    
    % Fourth pass - In RGCs with 3 inputs, see if we can re-assign 1 of
    % these inputs to an nearby RGC with only 1 cone input, ensuring that the cone types match
    minConeInputsNum = 3;
    maxInputsNum = 3;
    [rgcIDsWithThreeConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(connectionMatrix, minConeInputsNum, maxInputsNum, RGCRFPositionsMicrons);
    eccDegs =  WatsonRGCModel.rhoMMsToDegs(sqrt(sum(RGCRFPositionsMicrons(rgcIDsWithThreeConeInputs,:).^2,2))/1000.0);
    
    fprintf('\n -PASS 4: Reassigning 3-rd cone input of %d 3-cone input RGCs to a nearby orphanRGC or an RGC with 1 cone input ...',  numel(rgcIDsWithThreeConeInputs));
    tic
    
    successfullReassignementsNum = 0;
    orphanAssignmentsNum = 0;
    
    % Probabilities of cone-specific wiring
    pAttemptConeSpecificity = rand(1,numel(rgcIDsWithThreeConeInputs))<=maximizeConeSpecificity/100;
    
    for iRGC = 1:numel(rgcIDsWithThreeConeInputs)
        
        % Check whether to attempt cone-specific wiring in the center
        if (pAttemptConeSpecificity(iRGC) == 0)
            continue;
        end
        
        % Get the RGC index
        threeInputRGCindex = rgcIDsWithThreeConeInputs(iRGC);
        coneIndices = indicesOfConnectedCones{iRGC};
        
        % See if we can re-assign one of the 3 cone inputs to either an orphanRGC
        % or to a nearby RGC with 1 cone input, matching the cone type
        inputConeTypes = coneTypes(coneIndices);
        lConeIndices = find(inputConeTypes == LCONE_ID);
        mConeIndices = find(inputConeTypes == MCONE_ID);
                
        % Determine the order by which we will attempt cones
        if ((numel(lConeIndices) == 3) || (numel(mConeIndices) == 3))
            coneInputOrder = [1 2 3];
        else
            % mixed 3-cone input. Select the odd cone
            if (numel(lConeIndices) == 1)
                 coneInputOrder(1) = find(inputConeTypes == LCONE_ID);
            else
                 coneInputOrder(1) = find(inputConeTypes == MCONE_ID);
            end
            coneInputOrder = [coneInputOrder(1) setdiff([1 2 3], coneInputOrder(1))];
        end
        
      
        if (eccDegs(iRGC) <= 4.0)
            multiplier = 51;
        else
            multiplier = 7;
        end
        
        % First chance - seek to connect one cone input to a nearby orphan RGC
        % Since these are orphans look in a large region (41xlocal spacing)
        orphanRGCindex = findNeigboringOrphanRGC(threeInputRGCindex, ...
                    multiplier*RGCRFSpacingsMicrons(threeInputRGCindex), RGCRFPositionsMicrons, connectionMatrix);
                
        if (~isempty(orphanRGCindex))
            orphanAssignmentsNum = orphanAssignmentsNum + 1;
            indexOfConeToBeReassigned = coneIndices(coneInputOrder(1));

            % disconnect indexOfConeToBeReassigned from the mismatchedConeInputRGCindex
            % and connect it to the orphanRGCindex
            [connectionMatrix,  numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                indexOfConeToBeReassigned, threeInputRGCindex, ...
                indexOfConeToBeReassigned, orphanRGCindex,  ...
                connectionMatrix,  numberOfConeInputs); 
            % update the position of the previously orphan RGC
            RGCRFPositionsMicrons(orphanRGCindex,:) = conePositionsMicrons(indexOfConeToBeReassigned,:);
        else
            foundOneConeInputRGCwithMatchedConeIndex = false;
            % Go through all the cones and see which one can connect
            for kkk = 1:numel(coneInputOrder)
                if (foundOneConeInputRGCwithMatchedConeIndex)
                    continue;
                end
                coneInputNo = coneInputOrder(kkk);
                indexOfConeToBeReassigned = coneIndices(coneInputNo);
                typeOfConeToBeReassigned = coneTypes(indexOfConeToBeReassigned);

                matchedConeInputRGCindex = findNeigboringRGCWithOneConeInputOfMatchedType(...
                    indexOfConeToBeReassigned, typeOfConeToBeReassigned, ...
                    1.5*coneSpacingsMicrons(indexOfConeToBeReassigned), ...
                    coneTypes, conePositionsMicrons, connectionMatrix);

                if (~isempty(matchedConeInputRGCindex))
                    foundOneConeInputRGCwithMatchedConeIndex = true;
                    successfullReassignementsNum = successfullReassignementsNum + 1;
                    % disconnect indexOfConeToBeReassigned from the threeInputRGCindex
                    % and connect it to the matchedRGCindex
                    [connectionMatrix,  numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                        indexOfConeToBeReassigned, threeInputRGCindex, ...
                        indexOfConeToBeReassigned, matchedConeInputRGCindex, ...
                        connectionMatrix,  numberOfConeInputs); 
                end
            end % for kkk
        end
        
    end % iRGC
    
    
    fprintf('\nCompleted successfully for %d RGCs + %d orphan RGCs in %f minutes\n', successfullReassignementsNum, orphanAssignmentsNum, toc/60);
end


function [connectionMatrix, numberOfConeInputs] = ...
        performPass5(conePositionsMicrons, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs)
    
    % Final phase. For RGCs with several (>=4) cone inputs, see if we can
    % assign some of the most distant inputs to nearby RGCs with less cone inputs
    minConeInputsNum = 4;
    maxConeInputsNum = Inf;
    [rgcIDsWithSeveralConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(...
        connectionMatrix, minConeInputsNum, maxConeInputsNum, RGCRFPositionsMicrons);
    
    fprintf('\n \t -PASS 5: Reassigning farthest cone input of %d multi-input RGCs (>=4 cone inputs) to the closest RGCs with fewer cone inputs ...', numel(rgcIDsWithSeveralConeInputs));
    tic
        
    for phaseRepeat = 1:3
        
        if (phaseRepeat > 1)
            [rgcIDsWithSeveralConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(...
                connectionMatrix, minConeInputsNum, maxConeInputsNum, RGCRFPositionsMicrons);
        end
        
        conesReassignedInPhase(phaseRepeat) = 0;
        
        for iRGC = 1:numel(rgcIDsWithSeveralConeInputs)
            % Get the RGC index and its cone input indices
            multiInputRGCindex = rgcIDsWithSeveralConeInputs(iRGC);
            multiInputRGCconeIndices = indicesOfConnectedCones{iRGC};

            % Compute the centroid of all our cones
            centroidOfInputs = mean(conePositionsMicrons(multiInputRGCconeIndices,:),1);

            % Determine the furthest cone from the centroid
            dd = sqrt(sum((bsxfun(@minus, conePositionsMicrons(multiInputRGCconeIndices,:), centroidOfInputs)).^2,2));
            [~,idx] = max(dd);
            indexOfConeToBeReassigned = multiInputRGCconeIndices(idx(1));
            distanceOfConeToBeReassignedToParentMRGC = dd(idx(1));
            
            % Retrieve the neigboring RGCs to that cone and their cone inputs
            [neighboringRGCindices, neighboringRGCconeInputIndices] = indicesOfNeigbboringRGCsToTargetCone(...
                conePositionsMicrons(indexOfConeToBeReassigned,:), ...
                RGCRFPositionsMicrons, 1.2*RGCRFSpacingsMicrons(multiInputRGCindex), ...
                multiInputRGCindex, connectionMatrix);

            coneReassigned = false;
            for k = 1:numel(neighboringRGCindices)
                neighboringRGCindex = neighboringRGCindices(k);
                neighboringRGCconeInputIDs = neighboringRGCconeInputIndices{k};
                centroidOfInputsToNeighboringRGC = mean(conePositionsMicrons(neighboringRGCconeInputIDs,:),1);
                distanceOfConeToBeReassignedToNeighboringMRGC = sqrt(sum((conePositionsMicrons(indexOfConeToBeReassigned,:)-centroidOfInputsToNeighboringRGC).^2));
                
                if ((numel(neighboringRGCconeInputIDs)+1 <= numel(multiInputRGCconeIndices)) && ...
                        (coneReassigned == false) && ...
                        (distanceOfConeToBeReassignedToParentMRGC > distanceOfConeToBeReassignedToNeighboringMRGC) )
                    coneReassigned = true;
                    conesReassignedInPhase(phaseRepeat) = conesReassignedInPhase(phaseRepeat) + 1;
                    % disconnect indexOfConeToBeReassigned from the multiInputRGCindex
                    % and connect it to the neighboringRGCindex
                    [connectionMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                        indexOfConeToBeReassigned, multiInputRGCindex, ...
                        indexOfConeToBeReassigned, neighboringRGCindex, ...
                        connectionMatrix,  numberOfConeInputs); 
                end
            end % for k
        end % iRGC
    end % phase repeat
    
    fprintf('\nCompleted successfully for %d multi (>=4) -input RGCs in %f minutes\n', sum(conesReassignedInPhase), toc/60);
    %conesReassignedInPhase
end

function plotEccOfOrphanRGCs(numberOfConeInputs, RGCRFPositionsMicrons, subplotIndex)
    % Plot initial distribution of orphanRGCs
    orphanRGCindices = find(numberOfConeInputs == 0);
    orphanRGCeccentricitiesMicrons = sqrt(sum(RGCRFPositionsMicrons(orphanRGCindices,:).^2,2));
    orphanRGCeccentricitiesDegs = WatsonRGCModel.rhoMMsToDegs(orphanRGCeccentricitiesMicrons/1000.0);
    figure(100);
    subplot(5,1,subplotIndex)
    histogram(orphanRGCeccentricitiesDegs, 0:0.25:30)
    set(gca, 'XLim', [0 15], 'XTick', 0:1:15);
    xlabel('eccentricity (degs)');
    ylabel('orphan RGCs');
    title(sprintf('Pass %d', subplotIndex));
    drawnow;
end


function [neighboringRGCindices, neighboringRGCconeInputIndices] = findNeighboringRGCsAndTheirConeInputs(rgcIndex, ...
            localRGCRFSpacingMicrons, RGCRFPositionsMicrons, connectionMatrix)
        
    % Find nearby RGCs, which are potential re-allocators for this RGC's cone inputs
    neighboringRGCindices = indicesOfNeigboringRFs(rgcIndex, ...
        localRGCRFSpacingMicrons, ...
        RGCRFPositionsMicrons, ...
        'ascend');

    % Return cone indices
    neighboringRGCconeInputIndices = cell(1, numel(neighboringRGCindices));
    for k = 1:numel(neighboringRGCindices)
        theConeIndices = find(squeeze(connectionMatrix(:,neighboringRGCindices(k))) == 1);
        neighboringRGCconeInputIndices{k} = theConeIndices;
    end   
end

function orphanRGCindex = findNeigboringOrphanRGC(rgcIndex, localRGCSpacingMicrons, RGCRFPositionsMicrons, connectionMatrix)
    % Find this cone's neigboring RGCs
    neighboringRGCindices = indicesOfNeigboringRFs(...
        rgcIndex, ...
        localRGCSpacingMicrons, ...
        RGCRFPositionsMicrons, ...
        'ascend');
    
    % Go through each of the neighboring cones
    orphanRGCindex = [];
    foundNearbyOrhpahRGC = false;
    
    for kk = 1:numel(neighboringRGCindices)
        if (foundNearbyOrhpahRGC)
            continue;
        end
        rgcIndex = neighboringRGCindices(kk);
        coneInputIndices = find(squeeze(connectionMatrix(:,rgcIndex)) == 1);
        if (isempty(coneInputIndices))
            orphanRGCindex = rgcIndex;
            foundNearbyOrhpahRGC = true;
        end
    end % kk
end


function matchedRGCindex = findNeigboringRGCWithOneConeInputOfMatchedType(indexOfConeToBeReassigned, typeOfConeToBeReassigned, localConeSpacingMicrons, coneTypes, conePositionsMicrons, connectionMatrix)
    % Find this cone's neigboring cones
    neighboringConeIndices = indicesOfNeigboringRFs(...
        indexOfConeToBeReassigned, ...
        localConeSpacingMicrons, ...
        conePositionsMicrons, ...
        'ascend');

    % Go through each of the neighboring cones
    matchedRGCindex = [];
    foundRGCwithMatchedCone = false;
    
    for kk = 1:numel(neighboringConeIndices)
        if (foundRGCwithMatchedCone)
            continue;
        end
        coneIndex = neighboringConeIndices(kk);
        % Find the RGC that this cone is connected to
        candidateRGCindex = find(squeeze(connectionMatrix(coneIndex,:)) == 1);
        % See if this RGC has one only 1 cone and if so, if this cone
        % is matched to the typeOfConeToBeReassigned
        connectedConeIDs = find(squeeze(connectionMatrix(:,candidateRGCindex)) == 1);
        if (numel(connectedConeIDs) == 1) && (coneTypes(connectedConeIDs(1)) == typeOfConeToBeReassigned)
            % Got it
            foundRGCwithMatchedCone = true;
            matchedRGCindex = candidateRGCindex;
        end
    end % for kk
end

function [neighboringRGCindices, neighboringRGCconeInputIndices] = indicesOfNeigbboringRGCsToTargetCone(targetConePositionMicrons, rgcPositionsMicrons, localSpacingMicrons, sourceRGCindex, connectionMatrix)
    distances = sqrt(sum((bsxfun(@minus, rgcPositionsMicrons, targetConePositionMicrons).^2),2));
    neighboringRGCindices= find(distances <= localSpacingMicrons);
    [~,idx] = sort(distances(neighboringRGCindices), 'ascend');
    neighboringRGCindices = neighboringRGCindices(idx);
    % Exclude the source RGC index
    neighboringRGCindices = setdiff(neighboringRGCindices, sourceRGCindex);

    % Return cone indices
    neighboringRGCconeInputIndices = cell(1, numel(neighboringRGCindices));
    for k = 1:numel(neighboringRGCindices)
        theConeIndices = find(squeeze(connectionMatrix(:,neighboringRGCindices(k))) == 1);
        neighboringRGCconeInputIndices{k} = theConeIndices;
    end
    
end


function neighboringIndices = indicesOfNeigboringRFs(targetIndex, localSpacingMicrons, positionsMicrons, sortingOrder)
   distances = sqrt(sum((bsxfun(@minus, positionsMicrons, positionsMicrons(targetIndex,:)).^2),2));
   neighboringIndices = find(distances <= localSpacingMicrons);
   [d,idx] = sort(distances(neighboringIndices), sortingOrder);
   if ((d(1)~=0) && (strcmp(sortingOrder, 'ascend')))
       error('D is not zero: %f', d(1));
   else
       idx = idx(2:end);
   end
   if ((d(end) ~= 0) && (strcmp(sortingOrder, 'descend')))
       error('D is not zero: %f', d(end));
   else
       idx = idx(1:(end-1));
   end
   % Do not include the targer RF itself, which would have 0 distance
   neighboringIndices = neighboringIndices(idx);
end


function  [rgcIDsWithThreeOrMoreInputs, coneInputIDs] = findRGCsWithThreeOrMoreConeInputs(connectionMatrix, minNumberOfInputs, maxNumberOfInputs, RGCRFPositionsMicrons)
    
    coneInputsNumForAllRGCs = full(squeeze(sum(connectionMatrix,1)));
    
    % IDs of multi-input RGCs
    rgcIDsWithThreeOrMoreInputs = find(...
        (coneInputsNumForAllRGCs(:) >= minNumberOfInputs) & ...
        (coneInputsNumForAllRGCs(:) <= maxNumberOfInputs) );
    
    % Return indices sorted according to the RGC eccentricity
    ecc = sum(RGCRFPositionsMicrons(rgcIDsWithThreeOrMoreInputs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    rgcIDsWithThreeOrMoreInputs = rgcIDsWithThreeOrMoreInputs(idx);
    
    % Return the cone input IDs for each of these multi-input RGCs
    coneInputIDs = cell(1, numel(rgcIDsWithThreeOrMoreInputs));
    for k = 1:numel(rgcIDsWithThreeOrMoreInputs)
        % Get the indices of cones connected to this RGC
        rgcIndex = rgcIDsWithThreeOrMoreInputs(k);
        indicesOfConeInputs = find(squeeze(connectionMatrix(:, rgcIndex)) == 1);
        coneInputIDs{k} = indicesOfConeInputs;
    end
end


function [rgcIDs, coneInputIDs] = findRGCsWithTwoMismatchedConeInputs(connectionMatrix, RGCRFPositionsMicrons, coneTypes)
    
    rgcIDs = [];
    coneInputIDs = [];

    rgcIDsWithTwoInputs = find(squeeze(sum(connectionMatrix,1)) == 2);
    for k = 1:numel(rgcIDsWithTwoInputs)
        % Get the indices of cones connected to this RGC
        rgcIndex = rgcIDsWithTwoInputs(k);
        indicesOfConeInputs = find(squeeze(connectionMatrix(:, rgcIndex)) == 1);
        % add to list 
        if (coneTypes(indicesOfConeInputs(1)) ~= coneTypes(indicesOfConeInputs(2)) )
            rgcIDs = cat(2, rgcIDs, rgcIndex);
            coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
        end
    end
    if (isempty(coneInputIDs))
        return;
    end
    
    % Return indeces sorted according to the RGC eccentricity
    ecc = sum(RGCRFPositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~=2)
        error('Size must be 2');
    end
end

function  [connectionMatrix,  numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                    oldConeIndex, oldRGCindex, newConeIndex, newRGCindex,  ...
                    connectionMatrix, numberOfConeInputs)           

    % DISCONNECT oldConeIndex, from oldRGCindex
    connectionMatrix(oldConeIndex, oldRGCindex) = 0; % disconnect
    numberOfConeInputs(oldRGCindex) = numberOfConeInputs(oldRGCindex)-1;
    if (numberOfConeInputs(oldRGCindex) == 0)
        error('We removed all cone inputs from this RGC. This should never have haneppened');
    end
    
    % And CONNECT newConeIndex to newRGCindex
    connectionMatrix(newConeIndex, newRGCindex) = 1; % strength
    numberOfConeInputs(newRGCindex) = numberOfConeInputs(newRGCindex)+1;
end