function [connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
    connectConesToMidgetRGCRFcenters(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        orphanRGCpolicy, maximizeConeSpecificity, coneTypes, desiredConesToRGCratios, visualizeProcess)
    
    % Define constants
    LCONE_ID = 2;
    MCONE_ID = 3;
    SCONE_ID = 4;
    
    % Numbers of neurons
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    phaseString = sprintf('Connecting %d cones to %d mRGCRF centers: ', conesNum, rgcsNum);
    fprintf('\n\t %s-PASS 1: connecting each cone to their closest RGC ...', phaseString);
    tic
    
    % Number of cone inputs in the RF center of each mRGC
    numberOfConeInputs = zeros(1,rgcsNum);
    
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
       
        % Accumulate indices for sparse array construction 
        nonSconeIndices = cat(2, nonSconeIndices, iCone);
        closestRGCindices = cat(2, closestRGCindices, closestRGCIndex);
        distances = cat(2, distances, d);
        
        % Count number of cones in each RGC
        numberOfConeInputs(closestRGCIndex) = numberOfConeInputs(closestRGCIndex)+1;
    end % for iCone
    
    % Generate sparse connection and distance matrices
    connectionMatrix = sparse(nonSconeIndices, closestRGCindices, ones([1 numel(nonSconeIndices)]), conesNum, rgcsNum);  % connections
    distanceMatrix = sparse(nonSconeIndices, closestRGCindices, distances, conesNum, rgcsNum);                           % distances
    
    fprintf('Completed in %f minutes\n', toc/60);
    if (visualizeProcess)
        visualizeConnectivity(3, 'Pass 1', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios))
    end
    
 
    
    % Second pass
    % Some RGCs will have N > 2 cone inputs. See if we can assign inputs
    % 3 ... N to other RGCs in the neigborhood that have the minimal number of
    % inputs <= N-1. 
    
    % Find RGCs with more than 2 cone inputs
    rgcIDsWithMoreThanTwoInputs = find(squeeze(sum(connectionMatrix,1)) > 2);
    fprintf('\n \t %s-PASS 2: reassigning cones in %d RGCs (that have > 2 cone inputs) to neighboring RGCs ...', phaseString, numel(rgcIDsWithMoreThanTwoInputs));
    tic

    for iRGC = 1:numel(rgcIDsWithMoreThanTwoInputs)
        
        % Get the RGC index
        rgcIndex = rgcIDsWithMoreThanTwoInputs(iRGC);

        % Get the indices of cones connected to this RGC
        indicesOfMultipleConesConnectedToThisRGC = find(squeeze(connectionMatrix(:, rgcIndex)) == 1);

        % Sort the indices of connected cones according to their distance
        % from the RGC. We want to start re-assigning the longest distance
        % inputs first
        inputDistances = squeeze(distanceMatrix(indicesOfMultipleConesConnectedToThisRGC, rgcIndex));
        [~, idx] = sort(inputDistances, 'descend');
        indicesOfMultipleConesConnectedToThisRGC = indicesOfMultipleConesConnectedToThisRGC(idx);
        
        % Exclude the cone with the smallest distance (last one since we sorted in descending order). 
        % We dont want to assign this one to another RGC.
        indicesOfMultipleConesConnectedToThisRGC = indicesOfMultipleConesConnectedToThisRGC(1:end-1);
        
        % Go through each of these cones (and as long as we still have more than 2 connected cones)
        % and see if we can assign some to other RGCs
        localRGCSpacingMicrons = RGCRFSpacingsMicrons(rgcIndex);
        
        iCone = 0;
        while (numberOfConeInputs(rgcIndex) > 2) && (iCone < numel(indicesOfMultipleConesConnectedToThisRGC))
            iCone = iCone+ 1;
            % Retrieve cone index
            coneIndex = indicesOfMultipleConesConnectedToThisRGC(iCone);
            
            % Distance of this cone to current RGC
            distanceOFCurrentConeInputToParentRGC = distanceMatrix(coneIndex, rgcIndex);
            
            % Compute distance of this cone to all RGCs
            conePosMicrons = conePositionsMicrons(coneIndex,:);
            distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
            
            % Find indices of RGCs in neighborhood of this cone that also
            % are also not much more distant (up to 20% more) that this cone is to its parent RGC
            rgcIndicesWithinReach = find(...
                distances <= localRGCSpacingMicrons & ...
                distances <= 1.2*distanceOFCurrentConeInputToParentRGC ...
                );

            % Find which RGC to connect to. This will be the closest RGC that
            % has the minimal # of cones already connected to it, and which
            % is less that the current number of inputs to this RGC
            
            % The minimum number of cone connections for the RGCs within reach 
            neighboringRGCConnectionsNum = sum(squeeze(connectionMatrix(:, rgcIndicesWithinReach)),1);
            minConnections = min(neighboringRGCConnectionsNum);

            % Forget it if the inputs to the nearby RGCs are more numerous
            % than in the current RGC
            if (minConnections+1 >= numberOfConeInputs(rgcIndex))
                continue;
            end
            
            % Sort neigboring RGCs according to their distance to the cone
            [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');
            
            % Go through the neihboring RGCs from the closest to the
            % farthest and find the first one that has the min # of connections
            keepGoing = true; k = 0;
            reassignedRGCindex = nan;
            while (keepGoing) && (k < numel(sortedIndices))
                k = k + 1;
                theNeighboringRGCIndex = rgcIndicesWithinReach(sortedIndices(k));
                if (neighboringRGCConnectionsNum(sortedIndices(k)) == minConnections)
                    % That's the one, exit the while loop
                    keepGoing = false;
                    reassignedRGCindex = theNeighboringRGCIndex;
                end
            end

            if (~isnan(reassignedRGCindex))
                % Update the connection/ distance matrices
                newDistance = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons(reassignedRGCindex,:), conePosMicrons).^2),2));
                [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                    coneIndex, rgcIndex, ...                               % existing connection to be dissolved
                    coneIndex, reassignedRGCindex, newDistance, ... % new connection to be established
                    connectionMatrix, distanceMatrix, numberOfConeInputs); % old matrices
            end
        end % while
    end % for iRGC
    fprintf('Completed in %f minutes\n', toc/60);
    
    
    if (visualizeProcess)
        visualizeConnectivity(4, 'Pass 2', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios));
    end
    
    
    % Look for RGCs with exactly 2 cone inputs that are of different type and see if can re-assign/exchange
    % one these inputs to a neighboring RGC that has one cone input making sure we match the cone type
    if (maximizeConeSpecificity)
        
        fprintf('\n \t %s-PASS 3.1: finding cones in RGCs that receive a mixed 2 cone-input ...', phaseString);
        tic
        
        rgcIDsWithExactlyTwoInputs = find(squeeze(sum(connectionMatrix,1)) == 2);
        
        % Find indices of RGCs with 2 cone inputs of different type
        rgcIDsWithTwoConeInputsOfDifferentType = [];
        coneIDsToRGCswithMixedCenters = [];
        
        for iRGC = 1:numel(rgcIDsWithExactlyTwoInputs)
            % Get RGC index
            rgcIndex = rgcIDsWithExactlyTwoInputs(iRGC);

            % Get the indices of cones connected to this RGC
            indicesOfConesConnectedToThisRGC = find(squeeze(connectionMatrix(:, rgcIndex)) == 1);
            
            % Check to see if the input cones have mismatched types (i.e. L and M)
            if (coneTypes(indicesOfConesConnectedToThisRGC(1)) ~= coneTypes(indicesOfConesConnectedToThisRGC(2)))
                rgcIDsWithTwoConeInputsOfDifferentType = cat(2,rgcIDsWithTwoConeInputsOfDifferentType, rgcIndex);
                coneIDsToRGCswithMixedCenters = cat(2, coneIDsToRGCswithMixedCenters, indicesOfConesConnectedToThisRGC);
            end
        end
        fprintf('Completed in %f minutes\n', toc/60);
        
        fprintf('\n \t %s-PASS 3.2: reassigning cones in %d RGCs that receive a mixed 2 cone-input ...', phaseString, numel(rgcIDsWithTwoConeInputsOfDifferentType));
        tic
        
        % Try to reassign one of the 2 cones in the mixed 2-cone RGCs
        for iRGC = 1:numel(rgcIDsWithTwoConeInputsOfDifferentType)
            
            % Retrieve the 2 (mixed) cones that are connected to this RGC
            rgcIndex = rgcIDsWithTwoConeInputsOfDifferentType(iRGC);
            indicesOfConesConnectedToThisRGC = coneIDsToRGCswithMixedCenters(:,iRGC);
            coneIndexToBeDisconnected = indicesOfConesConnectedToThisRGC(1);
            
            % Find neighboring RGCs to the current RGC
            rgcPosMicrons =  RGCRFPositionsMicrons(rgcIndex,:);
            distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, rgcPosMicrons).^2),2));
            
            % See if there are nearby orphan RGCs we can reassign this cone
            % to. Since orphanRGCs will most likely be due to S-cones, look
            % in a large region, like:  5*RGCRFSpacingsMicrons(rgcIndex)
            rgcIndicesWithinReach = find(distances <= 2.0*RGCRFSpacingsMicrons(rgcIndex));
            nearbyOrphanRGCs = find(numberOfConeInputs(rgcIndicesWithinReach) == 0);
            
            if (~isempty(nearbyOrphanRGCs))
                candidateRGCindex = rgcIndicesWithinReach(nearbyOrphanRGCs(1));
                % Update the connection/ distance matrices
                newDistance = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons(candidateRGCindex,:), conePositionsMicrons(coneIndexToBeDisconnected,:)).^2),2));
                [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                        coneIndexToBeDisconnected, rgcIndex, ...                               % existing connection to be dissolved
                        coneIndexToBeDisconnected, candidateRGCindex, newDistance, ...         % new connection to be established
                        connectionMatrix, distanceMatrix, numberOfConeInputs);                 % old matrices
                %fprintf(2,'Assigned mixed cone input from one RGC to a nearby ORPHAN rgc !!\n');
            else
                % Ok, no nearby orphan RGC, lets see if we can reassign one
                % of these 2 cones to one of those
                foundCandidateRGCForTransferingMixedConeInput = false;
                for k = 1:2
                    % If we already assigned the first cone, skip the rest
                    if (foundCandidateRGCForTransferingMixedConeInput)
                        continue;
                    end
                    
                    coneIndexToBeDisconnected = indicesOfConesConnectedToThisRGC(k);
                    
                    % Find nearby RGCs to this cone
                    distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePositionsMicrons(coneIndexToBeDisconnected,:)).^2),2));
                    rgcIndicesWithinReach = find(distances <= 1.2*RGCRFSpacingsMicrons(rgcIndex));
                    
                    % Sort neigboring RGCs according to their distance to this RGC
                    [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');
                    candidateRGCindices = rgcIndicesWithinReach(sortedIndices);
                    
                    for cRGC = 1:numel(candidateRGCindices)
                        candidateRGCindex = candidateRGCindices(cRGC);
                        if ((candidateRGCindex == rgcIndex) || (foundCandidateRGCForTransferingMixedConeInput == true))
                            continue
                        end
                        if (numberOfConeInputs(candidateRGCindex) == 1)
                            % Determine the cone type of the candiate RGC
                            idx = find(squeeze(connectionMatrix(:, candidateRGCindex)) == 1);
                            if (coneTypes(idx) == coneTypes(coneIndexToBeDisconnected))
                                % Found it
                                foundCandidateRGCForTransferingMixedConeInput = true;
                                
                                % Update the connection/ distance matrices
                                newDistance = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons(candidateRGCindex,:), conePositionsMicrons(coneIndexToBeDisconnected,:)).^2),2));
                                [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                                    coneIndexToBeDisconnected, rgcIndex, ...                               % existing connection to be dissolved
                                    coneIndexToBeDisconnected, candidateRGCindex, newDistance, ...         % new connection to be established
                                    connectionMatrix, distanceMatrix, numberOfConeInputs);                 % old matrices
                               % fprintf('Assigned mixed cone input from one RGC to a nearby rgc !!\n');
                            end
                        end
                    end % for cRGC
                    
                end % k
            end 
        end % iRGC
    end % maximineConeSpecificity
    fprintf('Completed in %f minutes\n', toc/60);
        
    % Fourth pass - dealing with remaining orphan RGCs
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    
    fprintf('\n \t %s-PASS 4: Trying to reallocate cones to %d orphan RGCs ...', phaseString, numel(RGCswithZeroInputs));
    tic
    
    switch (orphanRGCpolicy) % valid options: {'remove', 'share input'}
        case 'remove'
            for k = 1:numel(RGCswithZeroInputs)
                fprintf('%s Removing orphan RGC #%d at position %2.1f,%2.1f\n', ...
                    phaseString, RGCswithZeroInputs(k), RGCRFPositionsMicrons(RGCswithZeroInputs(k),1), ...
                    RGCRFPositionsMicrons(RGCswithZeroInputs(k),2));
            end
            indicesToKeep = setdiff(1:size(RGCRFPositionsMicrons,1), RGCswithZeroInputs);
            RGCRFPositionsMicrons = RGCRFPositionsMicrons(indicesToKeep,:);
            RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(indicesToKeep);
            connectionMatrix = connectionMatrix(:,indicesToKeep);
            numberOfConeInputs = numberOfConeInputs(indicesToKeep);
            
        case 'share input'
            for k = 1:numel(RGCswithZeroInputs)
                rgcIndex = RGCswithZeroInputs(k);
                rgcPos = RGCRFPositionsMicrons(rgcIndex,:);
                % Find the closest cone that is not an S-cone
                validConeIndices = find(coneTypes < 4);
                distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons(validConeIndices,:), rgcPos).^2),2));
                [~, idx] = min(distances);
                iCone = validConeIndices(idx);
                connectionMatrix(iCone, rgcIndex) = 1;                            % strength
                distanceMatrix(iCone, rgcIndex) = distances(idx);             % distance
                numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)+1;
                fprintf('%s Provided shared cone input to orphan RGC #%d (from cone at position %2.1f, %2.1f)\n', ...
                    phaseString, rgcIndex, conePositionsMicrons(iCone,1), conePositionsMicrons(iCone,2));
            end
            
        case 'steal input'
            RGCswithRealZeroInputs = [];
            for iRGC = 1:numel(RGCswithZeroInputs)
                orphanRGCindex = RGCswithZeroInputs(iRGC);
                
                % Find neighboring RGCs
                rgcPosMicrons =  RGCRFPositionsMicrons(orphanRGCindex,:);
                distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, rgcPosMicrons).^2),2));
                rgcIndicesWithinReach = find(distances <= 1.4*RGCRFSpacingsMicrons(orphanRGCindex));
            
                % Sort neigboring RGCs according to their distance to this RGC
                [sortedDistances, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');

                % Exclude the RGC itself, which will have zero distance
                sortedIndices = sortedIndices(2:end);
                candidateRGCindices = rgcIndicesWithinReach(sortedIndices);
                
                % Determine the max cone inputs to the nearby RGCs
                coneInputsToNearbyRGCs = numberOfConeInputs(candidateRGCindices);
                maxConeInputsInTheNeigborhood = max(coneInputsToNearbyRGCs);
                
                if (maxConeInputsInTheNeigborhood > 2)
                    meanConeInputsInTheNeiborhood = floor(mean(coneInputsToNearbyRGCs));
                    % Go through the neirhboring cones and steal 1
                    % cone from each until the previously orphan RGC has meanConeInputsInTheNeiborhood
                    stolenCones = 0;
                    for cRGC = 1:numel(candidateRGCindices)
                        if (stolenCones >= meanConeInputsInTheNeiborhood)
                            continue;
                        end
                        candidateRGCindex = candidateRGCindices(cRGC);
                        if (numberOfConeInputs(candidateRGCindex) > 1)
                            % Find the cones connected to this RGC
                            indicesOfConesConnectedToTheCandidateRGC = find(squeeze(connectionMatrix(:, candidateRGCindex)) == 1);
                            % Find the cone closest to the orphanRGC location
                            distanceToCones = sqrt(sum((bsxfun(@minus, conePositionsMicrons(indicesOfConesConnectedToTheCandidateRGC,:), rgcPosMicrons).^2),2));
                            [~,coneIndexToSteal] = min(distanceToCones);

                            % Update the connection/ distance matrices
                            [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                                  coneIndexToSteal, candidateRGCindex, ...                               % existing connection to be dissolved
                                  coneIndexToSteal, orphanRGCindex, 0, ...         % new connection to be established
                                  connectionMatrix, distanceMatrix, numberOfConeInputs);                 % old matrices

                            stolenCones = stolenCones + 1;
                        end
                        
                    end
                    %fprintf(2,'Eliminated oprhan RGC by stealing %d cone inputs from nearby RGCs with > %d cone inputs!!\n', stolenCones, meanConeInputsInTheNeiborhood);
                    
                else % 2 or 1 max cone inputs in the neighborhood
                    switch (maxConeInputsInTheNeigborhood)
                        case 2
                            % Steel on input so we dont have an oprhan RGC here
                            foundCandidateRGCForTransferingMixedConeInput = false;
                            for cRGC = 1:numel(candidateRGCindices)
                                candidateRGCindex = candidateRGCindices(cRGC);
                                if ( (sum(connectionMatrix(:,candidateRGCindex),1) == maxConeInputsInTheNeigborhood) && ...
                                     (numberOfConeInputs(candidateRGCindex) > 1) && (foundCandidateRGCForTransferingMixedConeInput == false) )
                                     
                                    foundCandidateRGCForTransferingMixedConeInput = true;

                                     % Find the 2 cones connected to this RGC
                                     indexOfConesConnectedToTheCandidateRGC = find(squeeze(connectionMatrix(:, candidateRGCindex)) == 1);

                                     % Update the connection/ distance
                                     % matrices (by disconnecting the first
                                     % cone)
                                     coneIndexToSteal = indexOfConesConnectedToTheCandidateRGC(1);
                                     [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                                        coneIndexToSteal, candidateRGCindex, ...                    % existing connection to be dissolved
                                        coneIndexToSteal, orphanRGCindex, 0, ...                    % new connection to be established
                                        connectionMatrix, distanceMatrix, numberOfConeInputs);      % old matrices
                          
                                     %fprintf(2,'Eliminated oprhan RGC by stealing 1 cone input from a nearby RGC with 2 cone inputs!!\n');

                                end
                            end % cRGC
                        case 1
                            %fprintf('All %d RGCs in the neighbohood of this RGC had 1 cone input. RGC will remain an orphan.\n', numel(candidateRGCindices));
                            % All nearby RGCs have 1 input, so this has to remain
                            % an orphan RGC, so we will have to remove it
                            RGCswithRealZeroInputs = cat(2,RGCswithRealZeroInputs, orphanRGCindex);
                    end  % switch
                end
            end % iRGC
             
            % Remove any remaining orphan RGCs
            indicesToKeep = setdiff(1:size(RGCRFPositionsMicrons,1), RGCswithRealZeroInputs);
            RGCRFPositionsMicrons = RGCRFPositionsMicrons(indicesToKeep,:);
            RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(indicesToKeep);
            connectionMatrix = connectionMatrix(:,indicesToKeep);
            numberOfConeInputs = numberOfConeInputs(indicesToKeep);
            
        otherwise
            error('Unknown orpanRGCpolicy: ''%s''.', orphanRGCpolicy)
    end % Switch (ophanPolicy)
    
    
    fprintf('Completed in %2.1f minutes. There are still %d remaining orphan RGCs, which were removed.\n', toc/60, numel(RGCswithRealZeroInputs));         
    
    if (visualizeProcess)
        visualizeConnectivity(5, 'Pass 4', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios));
    end
   
    % Make sure we have no RGCs with 0 inputs
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    RGCswithZeroInputsFromConnectionMatrix = find(squeeze(sum(connectionMatrix,1)) == 0);
    fprintf('Number of RGCs with zero inputs: %d , from matrix:%d\n', numel(RGCswithZeroInputs), numel(RGCswithZeroInputsFromConnectionMatrix));
    
    for k = 1:numel(RGCswithZeroInputsFromConnectionMatrix)
        pos = RGCRFPositionsMicrons(RGCswithZeroInputs(k),:);
        error('RGC at position %2.1f, %2.1f has zero cone inputs\n', pos(1), pos(2));
    end
    
    if (visualizeProcess)
        visualizeConnectivity(10, 'Final pass', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix,  coneTypes, mean(desiredConesToRGCratios));
    end
end

function  [connectionMatrix, distanceMatrix, numberOfConeInputs] = updateConnectionAndDistanceMatrices(...
                    oldConeIndex, oldRGCindex, newConeIndex, newRGCindex, newDistance, ...
                    connectionMatrix, distanceMatrix, numberOfConeInputs)           

    % DISCONNECT oldConeIndex, from oldRGCindex
    connectionMatrix(oldConeIndex, oldRGCindex) = 0; % disconnect
    distanceMatrix(oldConeIndex, oldRGCindex) = 0; % distance
    numberOfConeInputs(oldRGCindex) = numberOfConeInputs(oldRGCindex)-1;
    if (numberOfConeInputs(oldRGCindex) == 0)
        error('We removed all cone inputs from RGC. This should never have haneppened');
    end
    
    % And CONNECT newConeIndex to newRGCindex
    connectionMatrix(newConeIndex, newRGCindex) = 1; % strength
    distanceMatrix(newConeIndex, newRGCindex) = newDistance;
    numberOfConeInputs(newRGCindex) = numberOfConeInputs(newRGCindex)+1;
end
