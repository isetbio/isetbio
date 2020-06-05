 % Phase 4: Connect cones to the mRGC RF centers
function runPhase4(runParams)

    % Load data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios');
        
    visualizeConnectionProcess = false;
    [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = computeConnectionMatrix(...
                RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
                coneTypes, desiredConesToRGCratios, runParams.orphanRGCpolicy, visualizeConnectionProcess);
            
    % Save connection data
    save(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.outputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios', 'midgetRGCconnectionMatrix');
end

function [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
    computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, ...
    RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
    coneTypes, desiredConesToRGCratios,  orphanRGCpolicy, visualizeProcess)
        
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
        orphanRGCpolicy, coneTypes, desiredConesToRGCratios, visualizeProcess);
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
    fprintf('%s Will align %d of %d RGCs in this patch, which had a cone-to-RGC ratio < 2\n', ...
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
    end
    
    fprintf('%s Could not find a cone that was not already aligned to another RGC for %d out of %d RGCs\n', phaseString, notAlignedRGCs, numel(indicesOfRGCsrequiringAlignment));
           
end

function [connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
    connectConesToMidgetRGCRFcenters(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        orphanRGCpolicy, coneTypes, desiredConesToRGCratios, visualizeProcess)
    
    phaseString = 'Connecting: ';
    
    % Numbers of neurons
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
     
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    maxNumberOfConnections = rgcsNum*5;
    
    fprintf('%s Allocating sparse arrays for connecting %d cones to %d RGCs\n', ...
        phaseString, conesNum, rgcsNum);

    connectionMatrix = spalloc(conesNum, rgcsNum, maxNumberOfConnections);
    distanceMatrix = spalloc(conesNum, rgcsNum, maxNumberOfConnections);
    numberOfConeInputs = zeros(1,rgcsNum);
    
    for iCone = 1:conesNum
        % Do not connect S-cones
        if (coneTypes(iCone) == 4)
            continue;
        end
        
        % Compute distance of this cone to all RGCs
        conePosMicrons = conePositionsMicrons(iCone,:);
        distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
        [~, connectedRGCIndex] = min(distances);
        
        connectionMatrix(iCone, connectedRGCIndex) = 1;                            % strength
        distanceMatrix(iCone, connectedRGCIndex) = distances(connectedRGCIndex);   % distance
        numberOfConeInputs(connectedRGCIndex) = numberOfConeInputs(connectedRGCIndex)+1;
    end % for iCone
    
    if (visualizeProcess)
        visualizeConnectivity(3, 'Pass 1', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios))
    end
    
    % Second pass
    % Some RGCs will have N > 2 cone inputs. See if we can assign inputs
    % 3 ... N to other RGCs in the neigborhood that have the minimal number of
    % inputs <= N-1
    % 
    rgcIDsWithMoreThanTwoInputs = find(squeeze(sum(connectionMatrix,1)) > 2);
    fprintf('%s There are %d out of a total of %d RGCs that receive more than 2 cone inputs.\n', ...
        phaseString, numel(rgcIDsWithMoreThanTwoInputs), size(connectionMatrix,2));

    for iRGC = 1:numel(rgcIDsWithMoreThanTwoInputs)
        
        % Get RGC index
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
        
        % Go throuch each of these cones and see if we can assign some
        % to other RGCs
        localRGCSpacingMicrons = RGCRFSpacingsMicrons(rgcIndex);
        
        iCone = 0;
        while (numberOfConeInputs(rgcIndex) > 2) && (iCone < numel(indicesOfMultipleConesConnectedToThisRGC))
            iCone = iCone+ 1;
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
            
            % Sort neigboring RGCs according to their distance to the cone
            [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');

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
                %fprintf('%s Disconnecting cone %d from rgc %d (which has %d inputs) and connecting it to rgc %d, which had %d inputs\n', ...
                %    phaseString, coneIndex, rgcIndex, numberOfConeInputs(rgcIndex), betterAlternativeRGCindex, numberOfConeInputs(betterAlternativeRGCindex));
                % Of DISCONNECT coneIndex from rgcIndex
                connectionMatrix(coneIndex, rgcIndex) = 0; % strength
                distanceMatrix(coneIndex, rgcIndex) = 0; % distance
                numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)-1;
                % And CONNECT coneIndex to betterAlternativeRGCindex
                connectionMatrix(coneIndex, betterAlternativeRGCindex) = 1; % strength
                distanceMatrix(coneIndex, betterAlternativeRGCindex) = distances(betterAlternativeRGCindex); % distance
                numberOfConeInputs(betterAlternativeRGCindex) = numberOfConeInputs(betterAlternativeRGCindex)+1;
            end
        end % while
    end % for iRGC
   
    if (visualizeProcess)
        visualizeConnectivity(4, 'Pass 2', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios));
    end
    
    
    % Third pass - dealing with orphan RGC (RGCs with no cone inputs)
    RGCswithZeroInputs = find(numberOfConeInputs == 0);
    
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
        otherwise
            error('Unknown orpanRGCpolicy: ''%s''.', orphanRGCpolicy)
    end
    
    if (visualizeProcess)
        visualizeConnectivity(5, 'Pass 3', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, mean(desiredConesToRGCratios));
    end
    

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
            coneIndicesConnectedToThisRGC = find(squeeze(connectionMatrix(:, theRGCindex)) == 1);
            
            % Adjust connection strengths
            distancesOfConesConnectedToThisRGC = squeeze(distancenMatrix(coneIndicesConnectedToThisRGC, theRGCindex));
            originalConnectionStrengths = squeeze(connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex));
            weights = exp(-0.2*(distancesOfConesConnectedToThisRGC/RGCRFSpacingsMicrons(theRGCindex)));
            adjustedConeInputWeights = originalConnectionStrengths .* weights;
            % Net weight of all inputs should be 1.
            adjustedConeInputWeights = adjustedConeInputWeights / sum(adjustedConeInputWeights);
            % Update connection matrix
            connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex) = adjustedConeInputWeights;
        end
    end
    
    if (visualizeProcess)
        visualizeConnectivity(10, 'Final pass', conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix,  coneTypes, mean(desiredConesToRGCratios));
    end
end