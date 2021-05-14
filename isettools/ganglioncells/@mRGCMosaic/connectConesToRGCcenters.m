% Method to connect cones within the FOV of the RGC mosaic to RGC RF centers
% Called by obj.wireRFcenterToInputCones()
function connectConesToRGCcenters(obj, coneRFPositionsMicrons, coneRFPositionsDegs, ...
    coneRFSpacingsMicrons, visualizeConnection)

    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass. Also
    % some RGCs will receive zero cone inputs because the closest cone was
    % an S-cone. This method sets the 
    % - obj.coneConnectivityMatrix and returns 
    % - the distance of each RGC to its closest cone.
    distances = connectEachConeToNearestRGC(obj,coneRFPositionsMicrons);
    
    if (visualizeConnection) && (1==2)
        figNo = 1;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after alignment');
    end
    
    % Report connectivity stats
    obj.connectivityStats(1);
    
    % Second pass. Minimize number of RGCs that connect to 3 cones by
    % recruiting orphan RGCs.
    [threeMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithThreeConeInputs(obj)  ; 
    minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
        coneRFPositionsMicrons, coneRFPositionsDegs, ...
        threeMismatchedInputConeTypeRGCs, coneInputIDs);
    
    if (visualizeConnection) && (1==2)
        figNo = 2;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of cones in 3-cone RGCS ');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(2);

    
    % Third pass. Minimize number of RGCs that connect to two mismatched cone
    % types by recruiting orphan RGCs
    
    % Find indices of all RGCs that have 2 mismatched cone inputs
    % These are returned in increasing eccentricity, so the first ones are
    % the most foveal ones. The coneInputIDs is a [N x 2] matrix with the 
    % indices of the 2 mismatced cone inputs to each of these RGCs
    [twoMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'mismatched');
    
    minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
        coneRFPositionsMicrons, coneRFPositionsDegs, ...
        twoMismatchedInputConeTypeRGCs, coneInputIDs);
    
    if (visualizeConnection) && (1==2)
        figNo = 3;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of cones in 2 mismatched-cone RGCS');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(3);
    
    
    % Fourth pass. Minimize number of RGCs that connect to two matched cone
    % types by recruiting orphan RGCs
    
    % Find indices of all RGCs that have 2 matched cone inputs
    % These are returned in increasing eccentricity, so the first ones are
    % the most foveal ones. The coneInputIDs is a [N x 2] matrix with the 
    % indices of the 2 mismatced cone inputs to each of these RGCs
    [twoMatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'matched');
    
    minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
        coneRFPositionsMicrons, coneRFPositionsDegs, ...
        twoMatchedInputConeTypeRGCs, coneInputIDs);
    
    if (visualizeConnection) && (1==2)
        figNo = 4;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of cones in 2 matched-cone RGCS');
    end
    
    % Report connectivity stats
    obj.connectivityStats(4);
    
end

function minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
    coneRFPositionsMicrons, coneRFPositionsDegs, multiInputRGCIndices, multiInputRGConeIndices)
    % Find RGCs with zero cone inputs. These are RGCs for which there was
    % no cone within the threshold distance, probably because the closest cone was an S-cone.
    orphanRGCIndices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 0);
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the RGC index
        mismatchedConeInputRGCindex = multiInputRGCIndices(iRGC);
        mismatchedConeIndices = multiInputRGConeIndices(iRGC,:);
            
        % Compute its radial eccentricity
        eccDegs2 = sum(obj.rgcRFpositionsDegs(mismatchedConeInputRGCindex,:).^2,2);
        
        % Search radius, equal to the cell's eccentricity with a min value of 0.5 degs
        minRadiusDegs = 0.5;
        searchRadiusDegs2 = max([minRadiusDegs^2 eccDegs2]);

        % Compute distances to all orphanRGCs 
        d2ToAllOrphanRGCs = sum((bsxfun(@minus, ...
            obj.rgcRFpositionsDegs(orphanRGCIndices,:), ...
            obj.rgcRFpositionsDegs(mismatchedConeInputRGCindex,:))).^2,2);
            
        [d2Min, idx] = min(d2ToAllOrphanRGCs);
        if (d2Min < searchRadiusDegs2)
            % OK, found orphanRGC within the search radius. We'll use it.
            theOrphanRGCindex = orphanRGCIndices(idx);
           
            % Remove it from the list of orphanRGCindices
            orphanRGCIndices = setdiff(orphanRGCIndices, theOrphanRGCindex);
            
            % Choose which cone input to reassign 
            if ((numel(mismatchedConeIndices) <= 2) || (mod(numel(mismatchedConeIndices),2) == 0))
                % Choose the first cone (arbitrary)
                indexOfConeToBeReassigned = mismatchedConeIndices(1);
            else
                % Choose the cone with the smallest population
                lconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.LCONE_ID);
                mconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.MCONE_ID);
                
                if (numel(lconeIndices) == 0)
                    fprintf('%d cone input RGC with %d lcones and %d mcones. Reassigning one M-cone\n', ...
                        numel(mismatchedConeIndices), numel(lconeIndices), numel(mconeIndices));
                    indexOfConeToBeReassigned = mismatchedConeIndices(mconeIndices(1));
                elseif (numel(mconeIndices) == 0)
                    fprintf('%d cone input RGC with %d lcones and %d mcones. Reassigning one L-cone\n', ...
                        numel(mismatchedConeIndices), numel(lconeIndices), numel(mconeIndices));
                    indexOfConeToBeReassigned = mismatchedConeIndices(lconeIndices(1));
                elseif (numel(lconeIndices) < numel(mconeIndices))
                    fprintf('%d cone input RGC with %d lcones and %d mcones. Reassigning one L-cone\n', ...
                        numel(mismatchedConeIndices), numel(lconeIndices), numel(mconeIndices));
                    indexOfConeToBeReassigned = mismatchedConeIndices(lconeIndices(1));
                else
                    fprintf('%d cone input RGC with %d lcones and %d mcones. Reassigning one M-cone\n', ...
                        numel(mismatchedConeIndices), numel(lconeIndices), numel(mconeIndices));
                    indexOfConeToBeReassigned = mismatchedConeIndices(mconeIndices(1));
                end
                
            end
            
            % Update the position of the previously orphan RGC to be that
            % of the cone that is reassigned to it
            obj.rgcRFpositionsMicrons(theOrphanRGCindex,:) = coneRFPositionsMicrons(indexOfConeToBeReassigned,:);
            obj.rgcRFpositionsDegs(theOrphanRGCindex,:) = coneRFPositionsDegs(indexOfConeToBeReassigned,:);

            % Note the obj.rgcRFspacingsMicrons, obj.rgcRFspacingsDegs
            % will now not reflect the local spacing 
            
            % Finally, update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  mismatchedConeInputRGCindex
            % and connecting 
            %   indexOfConeToBeReassigned  to theOrphanRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                indexOfConeToBeReassigned, mismatchedConeInputRGCindex, theOrphanRGCindex);
        else
            fprintf('Could not find an orphanRGC near 2-input RGC at location %2.2f,%2.2f\n', ...
                obj.rgcRFpositionsDegs(mismatchedConeInputRGCindex,1), ...
                obj.rgcRFpositionsDegs(mismatchedConeInputRGCindex,2));
        end
    end % iRGC
    
end

function updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, indexOfConeToBeReassigned, indexOfRGCToLooseCone, indexOfRGCtoReceiveCone)
    % DISCONNECT cone from its RGC
    if (obj.coneConnectivityMatrix(indexOfConeToBeReassigned, indexOfRGCToLooseCone) == 1)
        obj.coneConnectivityMatrix(indexOfConeToBeReassigned, indexOfRGCToLooseCone) = 0; % disconnect
    else
        error('Cone %d was not connected to RGC %d\n', indexOfConeToBeReassigned, indexOfRGCToLooseCone);
    end
    
    % And CONNECT it to the new RGC
    obj.coneConnectivityMatrix(indexOfConeToBeReassigned, indexOfRGCtoReceiveCone) = 1;
end


function [rgcIDs, coneInputIDs] = RGCsWithThreeConeInputs(obj)   
    % List of 3-cone input RGCs
    rgcIDs = [];
    % [N x 3] indices of the 3 cone input indices to each of the 3-cone RGC
    coneInputIDs = [];

    threeInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 3);
    for k = 1:numel(threeInputRGCindices)
        % Get the indices of cones connected to this RGC
        rgcIndex = threeInputRGCindices(k);
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        % Add to list 
        rgcIDs = cat(2, rgcIDs, rgcIndex);
        coneInputIDs = cat(2, coneInputIDs, indicesOfConeInputs);
    end
    if (isempty(coneInputIDs))
        return;
    end
    
    % Return 3-cone input RGC indices sorted according to the RGC eccentricity
    ecc = sum(obj.rgcRFpositionsMicrons(rgcIDs,:).^2,2);
    [~,idx] = sort(ecc, 'ascend');
    
    rgcIDs = rgcIDs(idx);
    coneInputIDs = coneInputIDs(:, idx);
    
    coneInputIDs = coneInputIDs';
    if (size(coneInputIDs,2) ~=3)
        error('Size must be 3');
    end
end


function [rgcIDs, coneInputIDs] = RGCsWithTwoConeInputs(obj, coneInputSchema)   
    % List of 2-cone input RGCs
    rgcIDs = [];
    % [N x 2] indices of the 2 cone input indices to each of the 2-cone RGC
    coneInputIDs = [];

    twoInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 2);
    for k = 1:numel(twoInputRGCindices)
        % Get the indices of cones connected to this RGC
        rgcIndex = twoInputRGCindices(k);
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


function distances = connectEachConeToNearestRGC(obj, coneRFPositionsMicrons)

    conesNum = size(coneRFPositionsMicrons,1);
    rgcsNum = size(obj.rgcRFpositionsMicrons,1);
    coneInputsNum = zeros(1,rgcsNum);
    
    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestRGCindices = [];
    nonSconeIndices = [];
    distances = [];
    
    for iCone = 1:conesNum
        % Do not connect S-cones to midget RGC centers
        if (obj.inputConeMosaic.coneTypes(iCone) == cMosaic.SCONE_ID)
            continue;
        end
        
        % Find the index of the closest RGC and connect the iCone to it
        [d2, nearestRGCIndex] = min(sum((bsxfun(@minus, obj.rgcRFpositionsMicrons, coneRFPositionsMicrons(iCone,:)).^2),2));
       
        d = sqrt(d2);
        maxDistanceBetweenConeAndRGC = 3.0*obj.rgcRFspacingsMicrons(nearestRGCIndex);
        if (d > maxDistanceBetweenConeAndRGC)
            %fprintf('Cone is too far from nearest RGC. Will not get connected to any RGC.\n');
            continue
        end
        
        % Accumulate indices for sparse array construction 
        nonSconeIndices = cat(2, nonSconeIndices, iCone);
        nearestRGCindices = cat(2, nearestRGCindices, nearestRGCIndex);
        
        % Distance of this RGC to its closest cone
        distances = cat(2, distances, d);
        
        % Update number of cones connected to this RGC
        coneInputsNum(nearestRGCIndex) = coneInputsNum(nearestRGCIndex)+1;
    end % iCone
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    obj.coneConnectivityMatrix = sparse(...
        nonSconeIndices, nearestRGCindices, ones([1 numel(nonSconeIndices)]), conesNum, rgcsNum);


%     for RGCindex = 1:rgcsNum    
%         connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, RGCindex)));
%         inputConeIDs = find(connectivityVector > 0.01);
%         if (numel(inputConeIDs) == 0)
%             orphanRGCIndices = cat(2, orphanRGCIndices, RGCindex);
%         end
%     end
    
end

