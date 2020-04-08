function [connectionMatrix, RGCRFPositionsMicrons] = connectConesToRGC(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        orphanRGCpolicy, desiredConesToRGCratios, visualizeProcess)
    
    % Numbers of neurons
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
     
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    connectionMatrix = zeros(conesNum, rgcsNum,2);
    numberOfConeInputs = zeros(1,rgcsNum);
    for iCone = 1:conesNum
        % Compute distance of this cone to all RGCs
        conePosMicrons = conePositionsMicrons(iCone,:);
        distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
        [~, connectedRGCIndex] = min(distances);
        
        connectionMatrix(iCone, connectedRGCIndex,1) = 1;                            % strength
        connectionMatrix(iCone, connectedRGCIndex,2) = distances(connectedRGCIndex); % distance
        numberOfConeInputs(connectedRGCIndex) = numberOfConeInputs(connectedRGCIndex)+1;
    end % for iCone
    

     if (visualizeProcess)
        visualizeConnectivity(3, 'Pass 1', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios))
     end
    
    
    
    % Second pass
    % Some RGCs will have N > 2 cone inputs. See if we can assign inputs
    % 3 ... N to other RGCs in the neigborhood that have the minimal number of
    % inputs <= N-1
    % 
    rgcIDsWithMoreThanTwoInputs = find(squeeze(sum(squeeze(connectionMatrix(:, :,1)),1)) > 2);
    fprintf('There are %d out of a total of %d RGCs that receive more than 2 cone inputs 1\n', ...
        numel(rgcIDsWithMoreThanTwoInputs), size(connectionMatrix,2));

    for iRGC = 1:numel(rgcIDsWithMoreThanTwoInputs)
        
        % Get RGC index
        rgcIndex = rgcIDsWithMoreThanTwoInputs(iRGC);

        % Get the indices of cones connected to this RGC
        indicesOfMultipleConesConnectedToThisRGC = find(squeeze(connectionMatrix(:, rgcIndex,1)) == 1);

        % Sort the indices of connected cones according to their distance
        % from the RGC. We want to start re-assigning the longest distance
        % inputs first
        inputDistances = squeeze(connectionMatrix(indicesOfMultipleConesConnectedToThisRGC, rgcIndex,2));
        [~, idx] = sort(inputDistances, 'descend');
        indicesOfMultipleConesConnectedToThisRGC = indicesOfMultipleConesConnectedToThisRGC(idx);
        
        % Exclude the cone with the smallest distance (last one since we sorted in descending order). 
        % We dont want to assign this one to another RGC.
        indicesOfMultipleConesConnectedToThisRGC = indicesOfMultipleConesConnectedToThisRGC(1:end-1);
        
        % Go throuch each of these cones and see if we can assign some
        % to other RGCs
        iCone = 0;
        while (numberOfConeInputs(rgcIndex) > 2) && (iCone < numel(indicesOfMultipleConesConnectedToThisRGC))
            iCone = iCone+ 1;
            coneIndex = indicesOfMultipleConesConnectedToThisRGC(iCone);
            
            % Distance of this cone to current RGC
            distanceOFCurrentConeInputToParentRGC = connectionMatrix(coneIndex, rgcIndex,2);
            
            % Compute distance of this cone to all RGCs
            conePosMicrons = conePositionsMicrons(coneIndex,:);
            distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
            
            % Find indices of RGCs in neighborhood of this cone that also
            % are less distant that this cone is to its parent RGC
            rgcIndicesWithinReach = find(...
                distances < 0.6*(max(RGCRFSpacingsMicrons)+coneSpacingsMicrons(coneIndex)) & ...
                distances <= 1.2*distanceOFCurrentConeInputToParentRGC ...
                );
            
            % Sort neigboring RGCs according to their distance to the cone
            [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');

            % Find which RGC to connect to. This will be the closest RGC that
            % has the minimal # of cones already connected to it, and which
            % is less that the current number of inputs to this RGC
            % The minimum number of cone connections for the RGCs within reach 
            neighboringRGCConnectionsNum = sum(squeeze(connectionMatrix(:, rgcIndicesWithinReach,1)),1);
            minConnections = min(neighboringRGCConnectionsNum);

            % Forget it if the inputs to the nearby RGCs are more numerous
            % than in the current RGC
            if (minConnections+1 >= numberOfConeInputs(rgcIndex))
                continue;
            end
            
            
            keepGoing = true; k = 0;
            betterAlternativeRGCindex = nan;
            while (keepGoing) && (k < numel(sortedIndices))
                k = k + 1;
                theNeighboringRGCIndex = rgcIndicesWithinReach(sortedIndices(k));
                if (neighboringRGCConnectionsNum(sortedIndices(k)) == minConnections)
                    keepGoing = false;
                    betterAlternativeRGCindex = theNeighboringRGCIndex;
                end
            end

            if (~isnan(betterAlternativeRGCindex))
                fprintf('Disconnecting cone %d from rgc %d (which has %d inputs) and connecting it to rgc %d, which had %d inputs\n', ...
                    coneIndex, rgcIndex, numberOfConeInputs(rgcIndex), betterAlternativeRGCindex, numberOfConeInputs(betterAlternativeRGCindex));
                % Of DISCONNECT coneIndex from rgcIndex
                connectionMatrix(coneIndex, rgcIndex, 1) = 0; % strength
                connectionMatrix(coneIndex, rgcIndex, 2) = 0; % distance
                numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)-1;
                % And CONNECT coneIndex to betterAlternativeRGCindex
                connectionMatrix(coneIndex, betterAlternativeRGCindex, 1) = 1; % strength
                connectionMatrix(coneIndex, betterAlternativeRGCindex, 2) = distances(betterAlternativeRGCindex); % distance
                numberOfConeInputs(betterAlternativeRGCindex) = numberOfConeInputs(betterAlternativeRGCindex)+1;
            end
        end % while
    end % for iRGC
   
    if (visualizeProcess)
        visualizeConnectivity(4, 'Pass 2', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios));
    end
    
    
    % Third pass - dealing with orphan RGC (RGCs with no cone inputs)
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    
    switch (orphanRGCpolicy) % valid options: {'remove', 'share input'}
        case 'remove'
            for k = 1:numel(RGCswithZeroInputs)
                fprintf('Removing orphan RGC #%d at position %2.1f,%2.1f\n', ...
                    RGCswithZeroInputs(k), RGCRFPositionsMicrons(RGCswithZeroInputs(k),1), ...
                    RGCRFPositionsMicrons(RGCswithZeroInputs(k),2));
            end
            indicesToKeep = setdiff(1:size(RGCRFPositionsMicrons,1), RGCswithZeroInputs);
            RGCRFPositionsMicrons = RGCRFPositionsMicrons(indicesToKeep,:);
            connectionMatrix = connectionMatrix(:,indicesToKeep,:);
            numberOfConeInputs = numberOfConeInputs(indicesToKeep);
            
        case 'share input'
            for k = 1:numel(RGCswithZeroInputs)
                rgcIndex = RGCswithZeroInputs(k);
                rgcPos = RGCRFPositionsMicrons(rgcIndex,:);
                % Find the closest cone
                distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPos).^2),2));
                [~, iCone] = min(distances);
                connectionMatrix(iCone, rgcIndex,1) = 1;                            % strength
                connectionMatrix(iCone, rgcIndex,2) = distances(iCone);             % distance
                numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)+1;
                fprintf('Provided shared cone input to orphan RGC #%d (from cone at position %2.1f, %2.1f)\n', ...
                    rgcIndex, conePositionsMicrons(iCone,1), conePositionsMicrons(iCone,2));
            end
        otherwise
            error('Unknown orpanRGCpolicy: ''%s''.', orphanRGCpolicy)
    end
    
    if (visualizeProcess)
        visualizeConnectivity(5, 'Pass 3', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios));
    end
    
    
    % Fourth pass. Some RGCs will have N = 2 cone inputs. If those inputs are from cones
    % of the same type we do nothing. If the inputs are from different cone
    % types see if we can assign the second cone to another RGC that has a
    % SINGLE input of the same cone type
    % NEED TO IMPLEMENT THIS
    

    % Make sure we have no RGCs with 0 inputs
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    for k = 1:numel(RGCswithZeroInputs)
        pos = RGCRFPositionsMicrons(RGCswithZeroInputs(k),:);
        error('RGC at position %2.1f, %2.1f has zero cone inputs\n', pos(1), pos(2));
    end
        
    normalizeNetWeights = ~true;
    if (normalizeNetWeights)
        % Some RGCs will have more than one cones connected to them. Adjust the
        % connection strengths of those cones depending on proximity and
        % cone-to-RGC ratio for each RGC
        % Connection matrix for RGCs with more than 1 cone inputs
        RGCswithMoreThanOneInputs = find(numberOfConeInputs > 1);
        for k = 1:numel(RGCswithMoreThanOneInputs)
            theRGCindex = RGCswithMoreThanOneInputs(k);
            % Find cones connected to this RGC
            coneIndicesConnectedToThisRGC = find(squeeze(connectionMatrix(:, theRGCindex, 1)) == 1);
            
            % Adjust connection strengths
            distancesOfConesConnectedToThisRGC = squeeze(connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex,2));
            originalConnectionStrengths = squeeze(connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex, 1));
            weights = exp(-0.2*(distancesOfConesConnectedToThisRGC/RGCRFSpacingsMicrons(theRGCindex)));
            adjustedConeInputWeights = originalConnectionStrengths .* weights;
            % Net weight of all inputs should be 1.
            adjustedConeInputWeights = adjustedConeInputWeights / sum(adjustedConeInputWeights);
            % Update connection matrix
            connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex, 1) = adjustedConeInputWeights;
        end
    end
    
    if (visualizeProcess)
        visualizeConnectivity(10, 'Final pass', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios));
    end
    
    % Return the connection strengths 
    connectionMatrix = squeeze(connectionMatrix(:,:,1));
end