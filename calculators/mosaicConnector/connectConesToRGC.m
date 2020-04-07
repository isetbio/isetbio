function connectionMatrix = connectConesToRGC(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess)
    
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
    
    % A small number of RGCs may have zero cone inputs as all local cones
    % have been assigned to some other RGCs. In this case select the
    % closest cone to provide shared input to those RGCs
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    for k = 1:numel(RGCswithZeroInputs)
        rgcIndex = RGCswithZeroInputs(k);
        rgcPos = RGCRFPositionsMicrons(rgcIndex,:);
        % Find the closest cone
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPos).^2),2));
        [~, iCone] = min(distances);
        connectionMatrix(iCone, rgcIndex,1) = 1;                            % strength
        connectionMatrix(iCone, rgcIndex,2) = distances(iCone);             % distance
        numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)+1;
        fprintf('Cone at position %2.1f, %2.1f provides shared input to more RGCs than RGC #%d.\n', ...
            conePositionsMicrons(iCone,1), conePositionsMicrons(iCone,2), rgcIndex);
    end
        
    
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

        % Exclude the cone with the smallest distance. We dont want to
        % assign this one to another RGC
        [~,closestConeIdx] = min(squeeze(connectionMatrix(indicesOfMultipleConesConnectedToThisRGC, rgcIndex,2)));
        indicesOfMultipleConesConnectedToThisRGC = setdiff(indicesOfMultipleConesConnectedToThisRGC, ...
            indicesOfMultipleConesConnectedToThisRGC(closestConeIdx));

        % Go throuch each of these cones and see if we can assign some
        % to other RGCs
        iCone = 0;
        while (numberOfConeInputs(rgcIndex) > 2) && (iCone < numel(indicesOfMultipleConesConnectedToThisRGC))
            iCone = iCone+ 1;
            coneIndex = indicesOfMultipleConesConnectedToThisRGC(iCone);
            % Compute distance of this cone to all RGCs
            conePosMicrons = conePositionsMicrons(coneIndex,:);
            distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
            % Find indices of RGCs in neighborhood of thiscone
            rgcIndicesWithinReach = find(distances < 0.6*(max(RGCRFSpacingsMicrons)+coneSpacingsMicrons(coneIndex)));
            % Sort neigboring RGCs according to their distance to the cone
            [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');

            % Find which RGC to connect to. This will be the closest RGC that
            % has the minimal # of cones already connected to it, and which
            % is less that the current number of inputs to this RGC
            % The minimum number of cone connections for the RGCs within reach 
            neighboringRGCConnectionsNum = sum(squeeze(connectionMatrix(:, rgcIndicesWithinReach,1)),1);
            minConnections = min(neighboringRGCConnectionsNum);

            if (minConnections >= numberOfConeInputs(rgcIndex))
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
    
    % Third pass. Some RGCs will have N = 2 cone inputs. If those inputs are from cones
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