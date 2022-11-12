% Method to connect cones within the FOV of the RGC mosaic to RGC RF centers
% Called by obj.wireRFcenterToInputCones()
function connectConesToRGCcenters(obj, idxConesInsideRGCmosaic, visualizeConnection)

    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass. Also
    % some RGCs will receive zero cone inputs because the closest cone was
    % an S-cone. This method sets the 
    % - obj.coneConnectivityMatrix and returns 
    % - the distance of each RGC to its closest cone.
    distances = connectEachConeToNearestRGC(obj,idxConesInsideRGCmosaic);
    
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
    
    % How far away to look for nearby RGCs for donating cones
    % This is times the mean cone spacing at that location
    searchRadiusFactor = 1.2;
    
    
    % Optimize central retina
    rangeDegs = [0 1];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    optimizeSubregion(obj, 1, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);

    rangeDegs = [1 2];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = true;
    optimizeSubregion(obj, 2, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);
    
    rangeDegs = [2 3];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = false;
    optimizeSubregion(obj, 3, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);

    
    rangeDegs = [3 5];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = false;
    optimizeSubregion(obj, 4, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);
    
    
    rangeDegs = [5 7];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = false;
    optimizeSubregion(obj, 5, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);
    
    rangeDegs = [7 Inf];
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs = false;
    optimizeSubregion(obj, 6, coneRFPositionsMicrons, coneRFPositionsDegs, ...
        coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs, ...
        useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection);
    
    
    hFig = figure(990);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[0 1 -0.4 0.4]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');

    hFig = figure(991);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[1 2 -0.4 0.4]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');
        
        
    hFig = figure(992);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[2 3 -0.5 0.5]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');
        
    hFig = figure(993);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[3 5 -0.6 0.6]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');
        
    hFig = figure(994);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[5 7 -0.7 0.7]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');
        
        
    hFig = figure(995);
    axesHandle = [];
    showConnectedCones = true;
    domain = 'microns';
    obj.visualizeConeMosaicTesselation(...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, ...
            'visualizationLimits', 0.8*[7 10 -1 1]*300, ...
            'figureHandle', hFig, ...
            'axesHandle', axesHandle, ...
            'plotTitle', 'final');
        
        
        
end

function optimizeSubregion(obj, subregionNo, coneRFPositionsMicrons, coneRFPositionsDegs, ...
    coneRFSpacingsMicrons, searchRadiusFactor, rangeDegs,  ...
    useOrphanRGCsToSplitMatchedConeMultiInputRGCs, visualizeConnection)
    
    fprintf(2, 'Optimizing region between %2.1f and %2.1f degrees\n', rangeDegs(1), rangeDegs(2));

    
    % ======================= DEALING WITH 3-INPUT RGCS =======================
    
    % Pass 1. Minimize number of RGCs that connect to 3 cones by
    % trying to donate one cone to a nearby RGC with one cone input only
    % The nearby cone type of the single input RGC must match the donated
    % cone type
    donatedConeMustMatch = true;
    [threeConeInputRGCs, coneInputIDs] = RGCsWithThreeConeInputs(obj, rangeDegs); 
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        coneRFSpacingsMicrons, searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC) with matched cone type reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    if (visualizeConnection) && (1==2)
        figNo = 12;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of 1 cone in 3-cone RGCS to nearby RGC with 1 cone input');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(1+subregionNo*10);

    
    % Pass 2. If we still have RGCs connecting to 3 cones, try to donate
    % one of them to a nearby orphan RGC
    [threeConeInputRGCs, coneInputIDs] = RGCsWithThreeConeInputs(obj, rangeDegs); 
    
    successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
        threeConeInputRGCs, coneInputIDs);
    fprintf('\t %d/%d of (3 cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    if (visualizeConnection) && (1==2)
        figNo = 2;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of 1 cone in 3-cone RGCS to oprhanRGCs');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(2+subregionNo*10);
    
    % Pass 3. Minimize number of RGCs that connect to 3 cones by
    % trying to donate one cone to a nearby RGC with one cone input only
    % The nearby cone type of the single input RGC does not have to match the donated cone type
    donatedConeMustMatch = false;
    [threeConeInputRGCs, coneInputIDs] = RGCsWithThreeConeInputs(obj, rangeDegs); 
    successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        coneRFSpacingsMicrons, searchRadiusFactor, threeConeInputRGCs, coneInputIDs, donatedConeMustMatch);
    
    fprintf('\t %d/%d (3 cone RGC -> 1 cone RGC) with un-matched cone type reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    % Report connectivity stats
    obj.connectivityStats(3+subregionNo*10);
    
    
    % ======================= DEALING WITH 2-INPUT RGCS =======================
    
    
    % Pass 4. Minimize number of RGCs that connect to two mismatched cone
    % types by re-assigning one of the cones to a nearby RGC with 1 cone input RGCs
    % of matched cone type
    
    % Find indices of all RGCs that have 2 mismatched cone inputs
    % These are returned in increasing eccentricity, so the first ones are
    % the most foveal ones. The coneInputIDs is a [N x 2] matrix with the 
    % indices of the 2 mismatced cone inputs to each of these RGCs
    [twoMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'mismatched', rangeDegs);
    
    successfullReallocations = minimizeFrequencyOf2MismatchedInputRGCSByDonatingToNearbyRGC(obj, ...
        coneRFSpacingsMicrons, searchRadiusFactor, ...
        twoMismatchedInputConeTypeRGCs, coneInputIDs);
    fprintf('\t %d/%d of (2 mismatched cone RGC -> nearby 1 cone RGC with matched cone) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    if (visualizeConnection) && (1==2)
        figNo = 4;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of cones in 2 mismatched-cone RGCS');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(4+subregionNo*10);
    
    
    % Pass 4. Minimize number of RGCs that connect to two mismatched cone
    % types by recruiting orphan RGCs
    
    % Find indices of all RGCs that have 2 mismatched cone inputs
    % These are returned in increasing eccentricity, so the first ones are
    % the most foveal ones. The coneInputIDs is a [N x 2] matrix with the 
    % indices of the 2 mismatced cone inputs to each of these RGCs
    [twoMismatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'mismatched', rangeDegs);
    
    successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
        twoMismatchedInputConeTypeRGCs, coneInputIDs);
    fprintf('\t %d/%d of (2 mismatched cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
    
    if (visualizeConnection) && (1==2)
        figNo = 5;
        axesHandle = [];
        showConnectedCones = true;
        domain = 'microns';
        obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
            coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
            obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
            showConnectedCones, domain, 'after re-assignment of cones in 2 mismatched-cone RGCS');
            
    end
    
    % Report connectivity stats
    obj.connectivityStats(5+subregionNo*10);
    
    if (useOrphanRGCsToSplitMatchedConeMultiInputRGCs)
        % Pass 5. Minimize number of RGCs that connect to two matched cone
        % types by recruiting orphan RGCs

        % Find indices of all RGCs that have 2 matched cone inputs
        % These are returned in increasing eccentricity, so the first ones are
        % the most foveal ones. The coneInputIDs is a [N x 2] matrix with the 
        % indices of the 2 mismatced cone inputs to each of these RGCs
        [twoMatchedInputConeTypeRGCs, coneInputIDs] = RGCsWithTwoConeInputs(obj, 'matched', rangeDegs);

        successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
            twoMatchedInputConeTypeRGCs, coneInputIDs);
        fprintf('\t %d/%d of (2 matched cone RGC -> orphan RGC) reassignments.\n', successfullReallocations(1), successfullReallocations(2));
        
        if (visualizeConnection) && (1==2)
            figNo = 6;
            axesHandle = [];
            showConnectedCones = true;
            domain = 'microns';
            obj.visualizeConeMosaicTesselation(figNo, axesHandle, ...
                coneRFPositionsMicrons, coneRFSpacingsMicrons, ...
                obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, ...
                showConnectedCones, domain, 'after re-assignment of cones in 2 matched-cone RGCS');
        end

        % Report connectivity stats
        obj.connectivityStats(6+subregionNo*10);
    end
    
end

function successfullReallocations = minimizeFrequencyOf2MismatchedInputRGCSByDonatingToNearbyRGC(obj, ...
        coneRFSpacingsMicrons, searchRadiusFactor, ...
        multiInputRGCIndices, multiInputRGConeIndices)
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the RGC index
        mismatchedConeInputRGCindex = multiInputRGCIndices(iRGC);
        mismatchedConeIndices = multiInputRGConeIndices(iRGC,:);
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(coneRFSpacingsMicrons(mismatchedConeIndices));
        
        % We have one L and one M cone.
        % Randomly select which one to try to re-assign first
        if (rand < 0.5)
            indexOfConeToBeReassigned = mismatchedConeIndices(1);
            indexOfOtherConeToBeReassigned = mismatchedConeIndices(2);
        else
            indexOfConeToBeReassigned = mismatchedConeIndices(2);
            indexOfOtherConeToBeReassigned = mismatchedConeIndices(1);
        end
        
        % Look for neignboring RGC this type of cone input
        maxConeInputsOfNearbyRGC = 1;
        targetConeType = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
        
        theTargetRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, targetConeType, ...
            indexOfConeToBeReassigned, searchRadiusMicrons);
        
        if (~isempty(theTargetRGCindex))
            % Update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  mismatchedConeInputRGCindex
            % and connecting 
            %   indexOfConeToBeReassigned  to theTargetRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                    indexOfConeToBeReassigned, mismatchedConeInputRGCindex, theTargetRGCindex);

            % Update successfullReallocations
            successfullReallocations = successfullReallocations + 1;
            continue;
        end
        
        % Look for neignboring RGC the other type of cone input
        maxConeInputsOfNearbyRGC = 1;
        targetConeType = obj.inputConeMosaic.coneTypes(indexOfOtherConeToBeReassigned);
        
        theTargetRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, targetConeType, ...
            indexOfOtherConeToBeReassigned, searchRadiusMicrons);
        
        if (~isempty(theTargetRGCindex))
            % Update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  mismatchedConeInputRGCindex
            % and connecting 
            %   indexOfConeToBeReassigned  to theTargetRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                    indexOfOtherConeToBeReassigned, mismatchedConeInputRGCindex, theTargetRGCindex);
            
            % Update successfullReallocations
            successfullReallocations = successfullReallocations + 1;
        end
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end


function successfullReallocations = minimizeFrequencyOf3InputRGCSByDonatingOneConeToNearbyRGC(obj, ...
        coneRFSpacingsMicrons, searchRadiusFactor, ...
        multiInputRGCIndices, multiInputRGConeIndices, donatedConeMustMatch)
    
    successfullReallocations = 0;
    
    if (isempty(multiInputRGCIndices))
        % Nothing to do
        successfullReallocations = [0 0];
        return;
    end
    
    for iRGC = 1:numel(multiInputRGCIndices)
        % Get the RGC index
        mismatchedConeInputRGCindex = multiInputRGCIndices(iRGC);
        mismatchedConeIndices = multiInputRGConeIndices(iRGC,:);
        
        % Find the non-matched cone type. This is the cone type we want
        % the nearbyRGC to have as its single input.
        [indexOfConeToBeReassigned, typeOfConeToBeReassigned] = indexOfConeWithSmallestPopulation(obj, ...
            donatedConeMustMatch, mismatchedConeIndices, mismatchedConeInputRGCindex);
        
        if (isempty(indexOfConeToBeReassigned))
            % This means that the cone selected would be at the RF center,
            % so dont reassign it as this would create a hole in the RF
            continue;
        end
        
        % If we do not care about matching the type of the donated cone, make its type []
        if (donatedConeMustMatch == false)
            typeOfConeToBeReassigned = [];
        end
        
        % Search radius
        searchRadiusMicrons = searchRadiusFactor*mean(coneRFSpacingsMicrons(mismatchedConeIndices));
        
        % Looking for a nearby RGC with one cone input only
        maxConeInputsOfNearbyRGC = 1;
        
        % Retrieve neighboring RGCs with no more than these cone inputs: maxConeInputsOfNearbyRGC,
        % with the majority of these inputs being of cone type: typeOfConeToBeReassigned
        theTargetRGCindex = indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
            maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
            indexOfConeToBeReassigned, searchRadiusMicrons);
            
        if (isempty(theTargetRGCindex))
            continue;
        end
        
        % Update the connectivityMatrix, by disconnecting
        %   indexOfConeToBeReassigned  FROM  mismatchedConeInputRGCindex
        % and connecting 
        %   indexOfConeToBeReassigned  to theTargetRGCindex
        updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                indexOfConeToBeReassigned, mismatchedConeInputRGCindex, theTargetRGCindex);
            
        % Update successfullReallocations
        successfullReallocations = successfullReallocations + 1;
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end


function successfullReallocations = minimizeFrequencyOfMultiInputRGCSByRecruitingNearbyOrphanRGCs(obj, ...
    multiInputRGCIndices, multiInputRGConeIndices)

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
            if (numel(mismatchedConeIndices) == 2)
                % Choose the first cone (arbitrary)
                indexOfConeToBeReassigned = mismatchedConeIndices(1);
            else
                donatedConeMustMatch = false;
                indexOfConeToBeReassigned = indexOfConeWithSmallestPopulation(obj, donatedConeMustMatch, mismatchedConeIndices, mismatchedConeInputRGCindex);
            end
            

            % Update the connectivityMatrix, by disconnecting
            %   indexOfConeToBeReassigned  FROM  mismatchedConeInputRGCindex
            % and connecting 
            %   indexOfConeToBeReassigned  to theOrphanRGCindex
            updateConnectivityMatrixByReassigningConeToDifferentRGC(obj, ...
                indexOfConeToBeReassigned, mismatchedConeInputRGCindex, theOrphanRGCindex);
            
            % Update successfullReallocations
            successfullReallocations = successfullReallocations + 1;
        
        end
    end % iRGC
    
    % Make it a percentage
    successfullReallocations = [successfullReallocations  numel(multiInputRGCIndices)];
end

function [indexOfConeToBeReassigned, typeOfConeToBeReassigned] = indexOfConeWithSmallestPopulation(obj, donatedConeMustMatch, mismatchedConeIndices, parentRGCindex)
    % Choose the cone with the smallest population, so if we have 2 Lcones
    % and 1 M cone, choose the M cone. But only if the 3 cones are
    % equidistant from the RGC RF center. If they are not, make sure we
    % do not pick the one that is at the RF center, as that would leave an
    % RF with a whole in its center

    if (donatedConeMustMatch)
        lconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.LCONE_ID);
        mconeIndices = find(obj.inputConeMosaic.coneTypes(mismatchedConeIndices) == cMosaic.MCONE_ID);
        if (numel(lconeIndices) == 0)
            indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(mconeIndices), parentRGCindex);
            typeOfConeToBeReassigned = cMosaic.MCONE_ID;
        elseif (numel(mconeIndices) == 0)
            indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(lconeIndices), parentRGCindex);
            typeOfConeToBeReassigned = cMosaic.LCONE_ID;
        elseif (numel(lconeIndices) < numel(mconeIndices))
            if (numel(lconeIndices) == 1)
                % See if this L-cone is closest to the RF center. If it is,
                % do not accept it
                targetConeDistanceFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(lconeIndices(1)),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                nonTargetConeDistancesFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(mconeIndices),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                if (targetConeDistanceFromRFcenter < 0.7*mean(nonTargetConeDistancesFromRFcenter))
                    % dont accept, too close to RF center
                    indexOfConeToBeReassigned = [];
                else
                    % accept it
                    indexOfConeToBeReassigned = mismatchedConeIndices(lconeIndices(1));
                end
            else
                indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(lconeIndices), parentRGCindex);
            end
            typeOfConeToBeReassigned = cMosaic.LCONE_ID;
        else
            if (numel(mconeIndices) == 1)
                % See if this M-cone is closest to the RF center. If it is,
                % do not accept it
                targetConeDistanceFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(mconeIndices(1)),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                nonTargetConeDistancesFromRFcenter = sum((bsxfun(@minus, ...
                    obj.inputConeMosaic.coneRFpositionsMicrons(mismatchedConeIndices(lconeIndices),:), ...
                    obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
                if (targetConeDistanceFromRFcenter < 0.7*mean(nonTargetConeDistancesFromRFcenter))
                    % dont accept, too close to RF center
                    indexOfConeToBeReassigned = [];
                else
                    % accept it
                    indexOfConeToBeReassigned = mismatchedConeIndices(mconeIndices(1));
                end
            else
                indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices(mconeIndices), parentRGCindex);
            end
            typeOfConeToBeReassigned = cMosaic.MCONE_ID;
        end
    else
        indexOfConeToBeReassigned = mostRemoteConeFromRFcenter(obj, mismatchedConeIndices, parentRGCindex);
        typeOfConeToBeReassigned = obj.inputConeMosaic.coneTypes(indexOfConeToBeReassigned);
    end
end

function coneIndex = mostRemoteConeFromRFcenter(obj, candidaceConeIndices, parentRGCindex)
    % Find the distances of these candidate cones from the RF center
    distancesFromRFcenter = sum((bsxfun(@minus, ...
            obj.inputConeMosaic.coneRFpositionsMicrons(candidaceConeIndices,:), ...
            obj.rgcRFpositionsMicrons(parentRGCindex,:))).^2,2);
    % Find cone that is most remote to the RF center
    [~,idx] = max(distancesFromRFcenter);
    coneIndex = candidaceConeIndices(idx);
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
    
    % Update the position of the donating RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, indexOfRGCToLooseCone)) == 1);
    obj.rgcRFpositionsMicrons(indexOfRGCToLooseCone,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(indexOfRGCToLooseCone,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
            
    % Update the position of the receiving RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, indexOfRGCtoReceiveCone)) == 1);
    obj.rgcRFpositionsMicrons(indexOfRGCtoReceiveCone,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(indexOfRGCtoReceiveCone,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
end


% Retrieve neighboring to target cone RGCs with no more than this many cone inputs: maxConeInputsOfNearbyRGC,
% with the majority of these inputs being of cone type: typeOfConeToBeReassigned
function targetRGCindex= indexOfNeighboringToConeRGCWithMaxConeInputs(obj, ...
       maxConeInputsOfNearbyRGC, typeOfConeToBeReassigned, ...
       indexOfConeToBeReassigned, searchRadiusMicrons)
   
    sourceConePositionMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(indexOfConeToBeReassigned,:);
    searchRadiusForNearbyRGCs = 1.5*searchRadiusMicrons;
    d2 = sum((bsxfun(@minus, obj.rgcRFpositionsMicrons, sourceConePositionMicrons)).^2,2);
    indicesOfRGCsWithinSearchRadius = find((d2 <= searchRadiusForNearbyRGCs^2)&(d2 > 0));
    
    if (numel(indicesOfRGCsWithinSearchRadius) > 1)
        % Sort them in ascending distance order
        [~,idx] = sort(d2(indicesOfRGCsWithinSearchRadius), 'ascend');
        indicesOfRGCsWithinSearchRadius = indicesOfRGCsWithinSearchRadius(idx);
    end
    

    targetRGCindex = [];
    % Go through all of these RGCs and measure the distances of their cone
    % inputs to the cone that is to be reassigned
    for iRGC = 1:numel(indicesOfRGCsWithinSearchRadius)
        
        % Get the rgc index
        rgcIndex = indicesOfRGCsWithinSearchRadius(iRGC);
        
        % Get the indices of cones connected to this RGC
        indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, rgcIndex)) == 1);
        
        % Determine if this RGC will meet the matching requirements 
        if (numel(indicesOfConeInputs) <= maxConeInputsOfNearbyRGC)
            
            % Find the closest of these cones to the source cone
            d2 = sum((bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:), sourceConePositionMicrons)).^2,2);
            [minDistance, idx] = min(d2);
            
            c1 = isempty(typeOfConeToBeReassigned)
            obj.inputConeMosaic.coneTypes(idx)
            
            c2 = (obj.inputConeMosaic.coneTypes(idx) == typeOfConeToBeReassigned)

            if (c1 || (~c1 && c2))
                % Add this to RGC and cone to the pool
                targetRGCindex = cat(2, targetRGCindex, rgcIndex);
                targetDistance = cat(2, targetDistance, minDistance);
            end
        end
    end
    
    % Pick the RGC with the min cone distance
    if (~isempty(targetRGCindex))
        idx = min(targetDistance);
        targetRGCindex = targetRGCindex(idx);
    end
    
end


        
function [rgcIDs, coneInputIDs] = RGCsWithThreeConeInputs(obj, rangeDegs)   
    % List of 3-cone input RGCs
    rgcIDs = [];
    % [N x 3] indices of the 3 cone input indices to each of the 3-cone RGC
    coneInputIDs = [];

    % Find all 3 input RGCs
    threeInputRGCindices = find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 3);
    
    % Select only those that lie between [minRadiusDegs, maxRadiusDegs]
    r2 = sum((obj.rgcRFpositionsDegs(threeInputRGCindices,:)).^2,2);
    idx = find((r2>=rangeDegs(1)^2) & (r2 <= rangeDegs(2)^2));
    threeInputRGCindices = threeInputRGCindices(idx);
    
    for k = 1:numel(threeInputRGCindices)
        rgcIndex = threeInputRGCindices(k);
        % Get the indices of cones connected to this RGC
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


function distances = connectEachConeToNearestRGC(obj, idxConesInsideRGCmosaic)


    % Note that these are the cones within the FOV of the mosaic
    % Not all the cones in the mosaic.
    conesNum = numel(idxConesInsideRGCmosaic);
    rgcsNum = size(obj.rgcRFpositionsMicrons,1);
    coneInputsNum = zeros(1,rgcsNum);
    
    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestRGCindices = [];
    nonSconeIndices = [];
    distances = [];
    
    for iCone = 1:conesNum
        % Cone index in the cMosaic
        coneIndex = idxConesInsideRGCmosaic(iCone);
        
        % Do not connect S-cones to midget RGC centers
        if (obj.inputConeMosaic.coneTypes(coneIndex) == cMosaic.SCONE_ID)
            continue;
        end
        
        % Find the index of the closest RGC and connect the cone to it
        [d2, nearestRGCIndex] = min(sum((bsxfun(@minus, obj.rgcRFpositionsMicrons, obj.inputConeMosaic.coneRFPositionsMicrons(coneIndex,:)).^2),2));
       
        d = sqrt(d2);
        maxDistanceBetweenConeAndRGC = 3.0*obj.rgcRFspacingsMicrons(nearestRGCIndex);
        if (d > maxDistanceBetweenConeAndRGC)
            %fprintf('Cone is too far from nearest RGC. Will not get connected to any RGC.\n');
            continue
        end
        
        % Accumulate indices for sparse array construction 
        nonSconeIndices = cat(2, nonSconeIndices, coneIndex);
        nearestRGCindices = cat(2, nearestRGCindices, nearestRGCIndex);
        
        % Distance of this RGC to its closest cone
        distances = cat(2, distances, d);
        
        % Update number of cones connected to this RGC
        coneInputsNum(nearestRGCIndex) = coneInputsNum(nearestRGCIndex)+1;
    end % iCone
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    conesNumFullMosaic = size(obj.inputConeMosaic.coneRFPositionsMicrons,1);
    obj.coneConnectivityMatrix = sparse(...
        nonSconeIndices, nearestRGCindices, ones([1 numel(nonSconeIndices)]), conesNumFullMosaic, rgcsNum);
    
end

