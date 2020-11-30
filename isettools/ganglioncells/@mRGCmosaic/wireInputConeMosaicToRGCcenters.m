function [coneConnectivityMatrix, rgcRFpositionsDegs, ...
    rgcRFpositionsMicrons, rgcRFspacingsDegs, rgcRFspacingsMicrons] = wireInputConeMosaicToRGCcenters(...
    rgcRFpositionsDegs, rgcRFpositionsMicrons,...
    conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    indicesOfConesNotConnectingToRGCcenters, coneSpecificityLevel, viewTesselationMaps)

            
    % Ensure that cones outside the fov of RGCmosaic do not get connected
    % to the RGC centers
    n = numel(indicesOfConesNotConnectingToRGCcenters);
    conePositionsDegs(indicesOfConesNotConnectingToRGCcenters,:) = Inf(n,2);
    conePositionsMicrons(indicesOfConesNotConnectingToRGCcenters,:) = Inf(n,2);
    coneSpacingsMicrons(indicesOfConesNotConnectingToRGCcenters) = Inf(n,1);
    
    % Compute local spacings from positions
    rgcRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(rgcRFpositionsMicrons);
    
    % Compute local densities from local spacings
    % Find the closest cone to each rgcRF
    [~,nearestConeIndices] = pdist2(conePositionsMicrons, rgcRFpositionsMicrons, 'euclidean', 'smallest',1);
    coneDensity = RGCmodels.Watson.convert.spacingToDensityForHexGrid(coneSpacingsMicrons(nearestConeIndices));
    rgcDensity = RGCmodels.Watson.convert.spacingToDensityForHexGrid(rgcRFspacingsMicrons);
    
    % Compute local ratio of cones to RGC from the local density ratios
    conesToRGCratios = coneDensity ./ rgcDensity;
    
    % Step1. Align each RGC with its nearest cone. This ensure all RGC's
    % are connected to at least one cone. Since cones are more numerous
    % than RGCs some cones will not connect to an RGC at this step. This 
    % step occurs only for RGCs for which the cone-to-RGC ratio is [1..2]
    visualizeProcess = ~true;
    [rgcRFpositionsMicrons, rgcRFpositionsDegs, rgcRFspacingsMicrons] = alignRGCmosaicToConeMosaic(...
         conePositionsDegs, ...
         conePositionsMicrons, ...
         coneSpacingsMicrons, ...
         coneTypes, ...
         rgcRFpositionsDegs, ...
         rgcRFpositionsMicrons, rgcRFspacingsMicrons, ...
         conesToRGCratios, visualizeProcess);
    
     % Step 2. Connect L and M cones to midget RGC centers.
     [coneConnectivityMatrix, rgcRFpositionsDegs, ...
      rgcRFpositionsMicrons, rgcRFspacingsDegs, rgcRFspacingsMicrons] = connectConesToMidgetRGCRFcenters(...
            conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            rgcRFpositionsDegs, rgcRFpositionsMicrons, rgcRFspacingsMicrons, ...
            coneSpecificityLevel, indicesOfConesNotConnectingToRGCcenters, ...
            viewTesselationMaps);

end


function [coneConnectivityMatrix, RGCRFPositionsDegs, RGCRFPositionsMicrons, ...
    RGCRFSpacingsDegs, RGCRFSpacingsMicrons] = ...
    connectConesToMidgetRGCRFcenters(conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes,...
        RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, coneSpecificityLevel, indicesOfConesNotConnectingToRGCcenters, viewTesselationMaps)
    
    
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    [coneConnectivityMatrix, numberOfConeInputs] = ...
        performPass1(conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, viewTesselationMaps);
    
    % Second pass - Maximize the frequency by which cones to 2-input RGCs have matched types
    [coneConnectivityMatrix, numberOfConeInputs, RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
        performPass2(conePositionsDegs, conePositionsMicrons,coneSpacingsMicrons, coneTypes, ...
        RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        coneConnectivityMatrix, numberOfConeInputs, coneSpecificityLevel, viewTesselationMaps);

    
    % Third pass. For RGCs with several (>=3) cone inputs, see if we can
    % assign some of the inputs to nearby RGCs with less cone inputs
    [coneConnectivityMatrix, numberOfConeInputs, RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
        performPass3(conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        coneConnectivityMatrix, numberOfConeInputs, viewTesselationMaps);


    % Fourth pass For RGCs with several (> 5) cone inputs, see if we can
    % assign some of the most distant inputs to nearby RGCs with less cone inputs
    [coneConnectivityMatrix, numberOfConeInputs] = ...
        performPass4(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        coneConnectivityMatrix, numberOfConeInputs, viewTesselationMaps);


    % Remove all remaining orphan RGCs
    orphanRGCindices = find(numberOfConeInputs == 0);
    orphanRGCindicesFromConnectionMatrix = find(squeeze(sum(coneConnectivityMatrix,1)) == 0);
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
    RGCRFPositionsDegs = RGCRFPositionsDegs(RGCindicesToKeep,:);
    RGCRFSpacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(RGCRFPositionsDegs);
    coneConnectivityMatrix = coneConnectivityMatrix(:,RGCindicesToKeep);
            

    fprintf('Checking connectivity matrix\n');
    assert(size(coneConnectivityMatrix,2) == size(RGCRFPositionsMicrons,1), ...
        'connectivity matrix size does not match size of positions');
    
    conesNum = size(conePositionsDegs,1);
    indicesOfConesConnectingToRGCcenters = setdiff(1:conesNum,indicesOfConesNotConnectingToRGCcenters);
    parfor iCone = 1:numel(indicesOfConesConnectingToRGCcenters)
        coneIndex = indicesOfConesConnectingToRGCcenters(iCone);
        checkConnectivityMatrix(coneConnectivityMatrix, coneTypes, conePositionsMicrons, coneIndex, 'at end')
    end
    
end

function checkConnectivityMatrix(connectionMatrix, coneTypes, conePositionsMicrons, targetConeIndex, prefix)
    % Check that targetCone (L- or M-) is connected to one and only one RGC
    if (coneTypes(targetConeIndex) ~= mRGCmosaic.SCONE_ID)
        rgcNumsConnectingTo(targetConeIndex) = full(sum(connectionMatrix(targetConeIndex,:),2));
        if (rgcNumsConnectingTo(targetConeIndex) ~= 1)
            fprintf('%s: cone at (%2.0f, %2.0f) is connected to %2.0f RGCs, NOT 1 RGC !! \n', ...
                prefix, conePositionsMicrons(targetConeIndex,1), conePositionsMicrons(targetConeIndex,2), rgcNumsConnectingTo(targetConeIndex));
        end
    end

end

function [connectionMatrix, numberOfConeInputs] = performPass1(...
    conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons, viewTesselationMaps)
    
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
        if (coneTypes(iCone) == mRGCmosaic.SCONE_ID)
            continue;
        end
        
        % Compute distance of this cone to all RGCs
        conePosMicrons = conePositionsMicrons(iCone,:);
        % Find the index of the closest RGC and connect the iCone to it
        [d, closestRGCIndex] = min(sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2)));
       
        maxDistanceBetweenConeAndRGC = 3.0*RGCRFSpacingsMicrons(closestRGCIndex);
        if (d > maxDistanceBetweenConeAndRGC)
            %fprintf('Cone is too far from nearest RGC. Will not get connected to any RGC.\n');
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

    % View tesselation
    if (viewTesselationMaps)
        figNo = 101;
        mRGCmosaic.renderTesselationMap(figNo, [], conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
    end
    fprintf('\nCompleted successfully in %f minutes\n',  toc/60);
end


    
function [connectionMatrix, numberOfConeInputs, RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
        performPass2(conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, coneSpecificityLevel, viewTesselationMaps)
    
    [rgcIDsWithTwoMismatchedConeInputs, indicesOfMismatchedCones] = findRGCsWithTwoMismatchedConeInputs(connectionMatrix, RGCRFPositionsMicrons,coneTypes);
    
    if (numel(rgcIDsWithTwoMismatchedConeInputs) == 0)
        fprintf('\n -PASS 2: Nothing to do. No RGCs with 2 mismatched cone inputs. \n');
    else
        
        if (viewTesselationMaps)
            % Clear figure
            figNo = 102;
            figure(figNo); clf;

             % View tesselation befor this step
            mRGCmosaic.renderTesselationMap(figNo, 1, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
        
        fprintf('\n -PASS 2: reassigning cones in %d RGCs with 2 mismatched cone inputs to neighboring RGCs with 1 or 0 cone inputs ...',  numel(rgcIDsWithTwoMismatchedConeInputs));
        tic
    
        successfullReassignementsNum = 0;
        orphanAssignmentsNum = 0;
    
        % Probabilities of cone-specific wiring
        pAttemptConeSpecificity = rand(1,numel(rgcIDsWithTwoMismatchedConeInputs))<=coneSpecificityLevel/100;
    
        for iRGC = 1:numel(rgcIDsWithTwoMismatchedConeInputs)
        
             % Check whether to attempt cone-specific wiring in the center
            if (pAttemptConeSpecificity(iRGC) == 0)
                continue;
            end
        
            % Get the RGC index
            mismatchedConeInputRGCindex = rgcIDsWithTwoMismatchedConeInputs(iRGC);
            mismatchedConeIndices = indicesOfMismatchedCones(iRGC,:);

            % Determine how far away to look for an orphan RGC
            rgcEccRadius = sqrt(sum(RGCRFPositionsDegs(mismatchedConeInputRGCindex,:).^2,2)); 
            if (rgcEccRadius <= 1.0)
                searchRadius = 71*RGCRFSpacingsMicrons(mismatchedConeInputRGCindex);
            else
                searchRadius = 7*RGCRFSpacingsMicrons(mismatchedConeInputRGCindex);
            end

            % Seek to connect one cone input to a nearby orphan RGC
            % Since these are orphans look in a large region (7xlocal spacing)
            orphanRGCindex = findNeigboringOrphanRGC(mismatchedConeInputRGCindex, ...
                searchRadius, RGCRFPositionsMicrons, connectionMatrix);

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
                RGCRFPositionsDegs(orphanRGCindex,:) = conePositionsDegs(indexOfConeToBeReassigned,:);
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
                        1.7*coneSpacingsMicrons(indexOfConeToBeReassigned), ...
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
    
    
    % Compute updated local spacings from updated positions
    RGCRFSpacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFPositionsMicrons);
    
    % View tesselation
    if (viewTesselationMaps)
        if (numel(rgcIDsWithTwoMismatchedConeInputs) > 0)
            mRGCmosaic.renderTesselationMap(figNo, 2, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
    end
end

function [connectionMatrix, numberOfConeInputs, RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = ...
        performPass3(conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        connectionMatrix, numberOfConeInputs, viewTesselationMaps)
    
    % Second pass. For RGCs with several (> 4) cone inputs, see if we can
    % assign some of the inputs to nearby RGCs with less cone inputs
    minConeInputsNum = 3;
    maxConeInputsNum = Inf;
    [rgcIDsWithSeveralConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(...
        connectionMatrix, minConeInputsNum, maxConeInputsNum, RGCRFPositionsMicrons);
    
    if (numel(rgcIDsWithSeveralConeInputs) == 0)
        fprintf('\n -PASS 3: Nothing to do (no RGCs with >= 3 cone inputs).\n');
    else
        
        if (viewTesselationMaps)
            % Clear figure
            figNo = 103;
            figure(figNo); clf;

            % View tesselation before this step
            mRGCmosaic.renderTesselationMap(figNo, 1, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
        
        fprintf('\n -PASS 3: Reassigning cone inputs of %d multi-input RGCs (>=3  cone inputs) to nearby RGCs with fewer cone inputs ...', numel(rgcIDsWithSeveralConeInputs));
        tic
    
        for iRGC = 1:numel(rgcIDsWithSeveralConeInputs)
            % Get the RGC index and its cone input indices
            multiInputRGCindex = rgcIDsWithSeveralConeInputs(iRGC);
            multiInputRGCconeIndices = indicesOfConnectedCones{iRGC};

            searchFactor = 1.3;
            % Retrieve all the neigboring RGCs and their cone input IDs
            [neighboringRGCindices, neighboringRGCconeInputIndices] = ...
                findNeighboringRGCsAndTheirConeInputs(multiInputRGCindex, ...
                searchFactor*RGCRFSpacingsMicrons(multiInputRGCindex), ...
                RGCRFPositionsMicrons, connectionMatrix);
            
            % Find the mean cone inputs in all neighbors
            meanInputsNumToNeighbors = 0;
            nonOrphanNeighbors = 0;
            for k = 1:numel(neighboringRGCindices)
                neighboringRGCconeInputIDs = neighboringRGCconeInputIndices{k};
                inputsNum = numel(neighboringRGCconeInputIDs);
                differenceInConeInputs(k) = numel(multiInputRGCconeIndices)-inputsNum;
                if (inputsNum>0)
                    meanInputsNumToNeighbors = meanInputsNumToNeighbors + inputsNum;
                    nonOrphanNeighbors = nonOrphanNeighbors + 1;
                end
            end
            meanInputsNumToNeighbors = meanInputsNumToNeighbors/nonOrphanNeighbors;
            
            % Scan each neighbor
            for k = 1:numel(neighboringRGCindices)
                neighboringRGCindex = neighboringRGCindices(k);
                neighboringRGCconeInputIDs = neighboringRGCconeInputIndices{k};

                if ((differenceInConeInputs(k)>1) && (meanInputsNumToNeighbors < numel(multiInputRGCconeIndices))) 
                    % Neighbor has less cone inputs and we still have more cone
                    % inputs than our neighbors on average

                    % Determine which cone input to give to this neighbor
                    if (numel(neighboringRGCconeInputIDs) == 0)
                        % Assigning to an orphanRGC, so just pick the first cone
                        indexOfConeToBeReassigned = multiInputRGCconeIndices(1);
                        % update the position of the previously orphan RGC
                        RGCRFPositionsMicrons(neighboringRGCindex,:) = conePositionsMicrons(indexOfConeToBeReassigned,:);
                        RGCRFPositionsDegs(neighboringRGCindex,:) = conePositionsDegs(indexOfConeToBeReassigned,:);
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
                    
                    % Update the neigbhor cone inputs
                    neighboringRGCconeInputIndices{k} = cat(1, neighboringRGCconeInputIndices{k}, indexOfConeToBeReassigned);
                end
            end % for k
        end % iRGC
        fprintf('\nCompleted successfully for %d multi (>=3) input RGCs in %f minutes\n', numel(rgcIDsWithSeveralConeInputs), toc/60);
    end
    
    % Compute updated local spacings from updated positions
    RGCRFSpacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFPositionsMicrons);
               
    if (viewTesselationMaps)
        if (numel(rgcIDsWithSeveralConeInputs) > 0)
            % View tesselation after this step
            mRGCmosaic.renderTesselationMap(figNo, 2, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
    end
    
end

function [connectionMatrix, numberOfConeInputs] = ...
        performPass4(conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
            connectionMatrix, numberOfConeInputs, viewTesselationMaps)
    
    % Final phase. For RGCs with several (>=4) cone inputs, see if we can
    % assign some of the most distant inputs to nearby RGCs with less cone inputs
    minConeInputsNum = 4;
    maxConeInputsNum = Inf;
    [rgcIDsWithSeveralConeInputs, indicesOfConnectedCones] = findRGCsWithThreeOrMoreConeInputs(...
        connectionMatrix, minConeInputsNum, maxConeInputsNum, RGCRFPositionsMicrons);
    
    if (numel(rgcIDsWithSeveralConeInputs) == 0)
         fprintf('\n -PASS 4: Nothing to do. No RGCs with 4+ cone inputs.\n');
    else
        
        if (viewTesselationMaps)
            % Clear figure
            figNo = 104;
            figure(figNo); clf;

            % View tesselation before this step
            mRGCmosaic.renderTesselationMap(figNo, 1, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                    RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
        
        fprintf('\n -PASS 4: Reassigning farthest cone input of %d multi-input RGCs (>=4 cone inputs) to the closest RGCs with fewer cone inputs ...', numel(rgcIDsWithSeveralConeInputs));
        tic

        conesReassignedInPhase = zeros(1,3);

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
    end
    
    if (viewTesselationMaps)
        % View tesselation after this step
        if (numel(rgcIDsWithSeveralConeInputs) > 0)
            mRGCmosaic.renderTesselationMap(figNo, 2, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                RGCRFPositionsMicrons, RGCRFSpacingsMicrons, connectionMatrix, 'microns');
        end
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


function  [RGCRFPositionsMicrons, RGCRFPositionsDegs, RGCRFSpacingsMicrons] = alignRGCmosaicToConeMosaic(...
    conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    RGCRFPositionsDegs, RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    desiredConesToRGCratios, visualizeProcess)
   
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
    
    % We only align mRGCs to L/M cones
    eligibleConeIndices = find(...
        (coneTypes == mRGCmosaic.LCONE_ID) | ...
        (coneTypes == mRGCmosaic.MCONE_ID) );
    
    conePositionsMicrons = conePositionsMicrons(eligibleConeIndices,:);
    conePositionsDegs = conePositionsDegs(eligibleConeIndices,:);
    
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
        rgcPDegs = RGCRFPositionsDegs(rgcIndex,:);
        rgcPMicrons = RGCRFPositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPMicrons).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 1.1*(RGCRFSpacingsMicrons(rgcIndex)));
        
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
        
        
        % If still all neighboring cones are already aligned to other RGCs, keep
        % original position
        if (isnan(alignmentConeIndex))
            notAlignedRGCs = notAlignedRGCs+1;
            %fprintf('%s Could not find a cone to align RGC %d that has not been aligned with another RGC already. Keeping original position %2.1f,%2.1f)\n', phaseString, rgcIndex, rgcPMicrons(1), rgcPMicrons(2));
            continue;
        end
        
        % Position of the cone to be aligned to
        rgcPMicronsConeAligned = conePositionsMicrons(alignmentConeIndex,:); 
        rgcPDegsConeAligned = conePositionsDegs(alignmentConeIndex,:); 
        
        % Update position of RGC to get closer to  rgcPMicronsConeAligned
        % depending on desiredConesToRGCratios(rgcIndex);
        g = min([1 desiredConesToRGCratios(rgcIndex)-1]);
        rgcPMicronsAligned = g * rgcPMicrons + (1-g)*rgcPMicronsConeAligned;
        rgcPDegsAligned = g * rgcPDegs + (1-g)*rgcPDegsConeAligned;
        
        % Update RGC mosaic
        RGCRFPositionMicronsBeforeAlignment = RGCRFPositionsMicrons(rgcIndex,:);
        RGCRFPositionsMicrons(rgcIndex,:) = rgcPMicronsAligned;
        RGCRFPositionsDegs(rgcIndex,:) = rgcPDegsAligned;
        
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
    
    % Update RGC RF spacings based on updated RGC RF positions
    RGCRFSpacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFPositionsMicrons);
    fprintf('%s Could not find a cone that was not already aligned to another RGC for %d out of %d RGCs\n', phaseString, notAlignedRGCs, numel(indicesOfRGCsrequiringAlignment));
           
end


% Support routines
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
    
    % Return indices sorted according to the RGC eccentricity
    ecc = sum(RGCRFPositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~=2)
        error('Size must be 2');
    end
end


function [neighboringRGCindicesSorted, neighboringRGCconeInputIndicesSorted] = findNeighboringRGCsAndTheirConeInputs(rgcIndex, ...
            localRGCRFSpacingMicrons, RGCRFPositionsMicrons, connectionMatrix)
        
    % Find nearby RGCs, which are potential re-allocators for this RGC's cone inputs
    neighboringRGCindices = indicesOfNeigboringRFs(rgcIndex, ...
        localRGCRFSpacingMicrons, ...
        RGCRFPositionsMicrons, ...
        'ascend');

    % Return cone indices
    neighboringRGCconeInputIndices = cell(1, numel(neighboringRGCindices));
    neighboringRGCconeInputsNum = zeros(1, numel(neighboringRGCindices));
    for k = 1:numel(neighboringRGCindices)
        theConeIndices = find(squeeze(connectionMatrix(:,neighboringRGCindices(k))) == 1);
        neighboringRGCconeInputIndices{k} = theConeIndices;
        neighboringRGCconeInputsNum(k) = numel(theConeIndices);
    end   
    
    % Sort neighboring RGCs according to the number of cone inputs
    [~,iidx] = sort(neighboringRGCconeInputsNum, 'ascend');
    neighboringRGCindicesSorted = neighboringRGCindices(iidx);
    neighboringRGCconeInputIndicesSorted = cell(1, numel(neighboringRGCindices));
    for k = 1:numel(iidx)
        neighboringRGCconeInputIndicesSorted{k} = neighboringRGCconeInputIndices{iidx(k)};
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


function  [rgcIDsWithThreeOrMoreInputs, coneInputIDs] = findRGCsWithThreeOrMoreConeInputs(connectionMatrix, ...
    minNumberOfInputs, maxNumberOfInputs, RGCRFPositionsMicrons)
    
    coneInputsNumForAllRGCs = full(squeeze(sum(connectionMatrix,1)));
    
    % IDs of multi-input RGCs
    rgcIDsWithThreeOrMoreInputs = find(...
        (coneInputsNumForAllRGCs(:) >= minNumberOfInputs) & ...
        (coneInputsNumForAllRGCs(:) <= maxNumberOfInputs) );
    
    sortPriority = 'eccentricity';
    sortPriority = 'number of inputs';
    
    switch (sortPriority)
        case 'eccentricity'
            % Return indices sorted according to the RGC eccentricity
            ecc = sum(RGCRFPositionsMicrons(rgcIDsWithThreeOrMoreInputs,:).^2,2);
            [~,idx] = sort(ecc, 'ascend');
            rgcIDsWithThreeOrMoreInputs = rgcIDsWithThreeOrMoreInputs(idx);
        case 'number of inputs'
            % Return indices sorted according to the number of cone inputs,
            % RGCs with most inputs first
            [~,idx] = sort(coneInputsNumForAllRGCs(rgcIDsWithThreeOrMoreInputs), 'descend');
            rgcIDsWithThreeOrMoreInputs = rgcIDsWithThreeOrMoreInputs(idx);
        otherwise
            error('Unknown sort priority: ''%s''.', sortPriority);
    end
    % Return the cone input IDs for each of these multi-input RGCs
    coneInputIDs = cell(1, numel(rgcIDsWithThreeOrMoreInputs));
    for k = 1:numel(rgcIDsWithThreeOrMoreInputs)
        % Get the indices of cones connected to this RGC
        rgcIndex = rgcIDsWithThreeOrMoreInputs(k);
        indicesOfConeInputs = find(squeeze(connectionMatrix(:, rgcIndex)) == 1);
        coneInputIDs{k} = indicesOfConeInputs;
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
