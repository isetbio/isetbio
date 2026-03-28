function theOptimizedDestinationRFindices = transferSourceRFsBetweenNearbyDestinationRFs(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('swapSourceRFsInsteadOffUnilateralTransfer', false, @islogical);
    p.addParameter('figureDir', fullfile(isetbioRootPath,'local'), @ischar);

    p.parse(varargin{:});
    swapSourceRFsInsteadOffUnilateralTransfer = p.Results.swapSourceRFsInsteadOffUnilateralTransfer;
    figureDir = p.Results.figureDir;
    generateProgressVideo = p.Results.generateProgressVideo;

    if (swapSourceRFsInsteadOffUnilateralTransfer)
        maxConesNumInTargetedRGCs = obj.wiringParams.maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs;
        convergenceSequenceFigNo = 5051;
    else
        maxConesNumInTargetedRGCs = obj.wiringParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs;
        convergenceSequenceFigNo = 5050;
    end

    debugConnectivityMatrix = false;

    % Compute the current costs for all destination RF to maintain their current inputs 
    theCostComponentsMatrix = obj.totalPoolingCosts();
    netCostSequences(1,:) = theCostComponentsMatrix;

    netTransfers = []; theOptimizedDestinationRFindices = [];
    for targetSourceRFsNum = 1:maxConesNumInTargetedRGCs
        if (swapSourceRFsInsteadOffUnilateralTransfer)
            fprintf('Optimizing swapping of source RFs in %d-input destination RFs...\n', targetSourceRFsNum);
            minNumberOfConeInputsInDonorNearbyRGCs = [];
            pdfFigureName = fullfile(figureDir, 'swapOptimization.pdf');
            phaseDescriptor = sprintf('swapping sourceRFs in %d input RFs',targetSourceRFsNum);
        else
            fprintf('Optimizing transfer of source RFs to %d-input destination RFs...\n', targetSourceRFsNum);
            minNumberOfConeInputsInDonorNearbyRGCs = 2 + max([targetSourceRFsNum maxConesNumInTargetedRGCs]);
            pdfFigureName = fullfile(figureDir, 'transferOptimization.pdf');
            phaseDescriptor = sprintf('transfering sourceRFs to %d input RFs', targetSourceRFsNum);
        end

        % Do a multi-pass
        [optimizedTargetRGCindices, netCostSequences, netTransfers, phaseDescriptor] = multiPassLocalConeInputAnomalyReduction(obj, targetSourceRFsNum, ...
            minNumberOfConeInputsInDonorNearbyRGCs, netCostSequences, netTransfers, convergenceSequenceFigNo, phaseDescriptor, debugConnectivityMatrix);

        theOptimizedDestinationRFindices = cat(1, theOptimizedDestinationRFindices(:), optimizedTargetRGCindices(:));
    end
    theOptimizedDestinationRFindices = unique(theOptimizedDestinationRFindices);

    hFig = figure(convergenceSequenceFigNo);
    NicePlot.exportFigToPDF(pdfFigureName, hFig, 300);

    % Save the metaDataStuct for this stage
    if (obj.saveIntermediateConnectivityStagesMetaData)
        obj.updateIntermediateMetaDataStructs(phaseDescriptor, netCostSequences, netTransfers);
    end
    
    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.intermediateFigureHandles{numel(obj.intermediateFigureHandles)+1} = ...
            obj.visualizeCurrentConnectivity(1003);
    end
end

function [allOptimizedDestinationRFindices, netCostSequences, netTransfers, phaseDescriptor] = multiPassLocalConeInputAnomalyReduction(obj, ...
    targetSourceRFsNum, minNumberOfConeInputsInDonorNearbyRGCs, netCostSequences, netTransfers, convergenceSequenceFigNo, phaseDescriptor, debugConnectivityMatrix)

    if (~isempty(minNumberOfConeInputsInDonorNearbyRGCs))
        theExaminedMinNumberOfSourceInputs = minNumberOfConeInputsInDonorNearbyRGCs:-1:(targetSourceRFsNum+2);
    else
        theExaminedMinNumberOfSourceInputs = 0;
    end

    allOptimizedDestinationRFindices = [];

    for minNumberOfSourceInputsToConsideredNearbyDestinationRFs = theExaminedMinNumberOfSourceInputs

        if (minNumberOfSourceInputsToConsideredNearbyDestinationRFs > 0)
            fprintf('\tExamining candidate transfers from %d+ input nearby destination RFs ...\n', ...
                minNumberOfSourceInputsToConsideredNearbyDestinationRFs);
            phaseDescriptor = sprintf('%s: from %d input RFs', phaseDescriptor, minNumberOfSourceInputsToConsideredNearbyDestinationRFs);
        end

        redoOptimization = true; currentPass = 0;
        while (redoOptimization) && (currentPass < obj.wiringParams.maxPassesNum)
            % Update currentPass
            currentPass = currentPass+1;

            % Do a single-pass anomaly reduction
            [destinationRFsNumWithTargetSourceRFsNum, optimizationTransfersNum, ...
             theOptimizedTargetDestinationRFindices] = singlePassLocalConeInputAnomalyReduction(obj, ...
                targetSourceRFsNum, minNumberOfSourceInputsToConsideredNearbyDestinationRFs, currentPass, debugConnectivityMatrix);

            % Accumulate theOptimizedTargetDestinationRFindices
            allOptimizedDestinationRFindices = cat(1, allOptimizedDestinationRFindices, theOptimizedTargetDestinationRFindices(:));

            if (optimizationTransfersNum == 0) || (destinationRFsNumWithTargetSourceRFsNum == 0)
                redoOptimization = false;
            else
                if (minNumberOfSourceInputsToConsideredNearbyDestinationRFs > 0)
                    fprintf('\t%d optimized transfer(s) in %d-input RFs to nearby %d+ input RFs during pass #%d (max passes: %d)\n', ...
                    optimizationTransfersNum, minNumberOfSourceInputsToConsideredNearbyDestinationRFs, targetSourceRFsNum, ...
                    currentPass, obj.wiringParams.maxPassesNum);
                    plotTitle = sprintf('source RFs num (transfer): %d -> %d',minNumberOfSourceInputsToConsideredNearbyDestinationRFs, targetSourceRFsNum);
                else
                    fprintf('%d optimized swap(s) in %d-input destination RFs with nearby destination RFs during pass #%d (max passes: %d)\n', ...
                    optimizationTransfersNum, targetSourceRFsNum, ...
                    currentPass, obj.wiringParams.maxPassesNum);
                    plotTitle = sprintf('source RFs num (swap): %d',targetSourceRFsNum);
                end
            end

            if (optimizationTransfersNum > 0)
                % Update transfers sequence
                netTransfers(numel(netTransfers)+1) = optimizationTransfersNum;

                % Compute the current costs for all destination RF to maintain their current inputs 
                theCostComponentsMatrix = obj.totalPoolingCosts();
                netCostSequences(size(netCostSequences,1)+1,:) = theCostComponentsMatrix;

                % Visualize convergence
                MosaicConnector.visualizeConvergenceSequence(currentPass, ...
                    netCostSequences, obj.costComponentNames(), ...
                    netTransfers, obj.wiringParams.maxPassesNum, ...
                    plotTitle, ...
                    convergenceSequenceFigNo);
            end %  if (optimizationTransfersNum > 0)

        end % while loop
    end % for minNumberOfSourceInputsToConsideredNearbyDestinationRFs 
end

function checkConnectivityMatrix(connectivityMatrix, targetSourceRFsNum, currentPass)

    destinationRFconvergenceVector = full(sum(connectivityMatrix,1));
    sourceDivergenceVector = full(sum(connectivityMatrix,2));
    fprintf('Optimizing RGCs with %d inputs. Pass %d.\n', targetSourceRFsNum, currentPass);
    fprintf('destinationRF convergence range: %d-%d\n', min(destinationRFconvergenceVector(:)), max(destinationRFconvergenceVector(:)));
    fprintf('sourceRF divergence range: %d-%d\n', min(sourceDivergenceVector(:)), max(sourceDivergenceVector(:)));

    if (max(sourceDivergenceVector(:))>1)
        error('cones connecting to more than 1 RGC')
    end
end

function [destinationRFsNumWithTargetSourceRFsNum, optimizationTransfersNum, theOptimizedTargetDestinationRFindices] = ...
        singlePassLocalConeInputAnomalyReduction(obj, targetSourceRFsNum, minNumberOfSourceInputsToConsideredNearbyDestinationRFs,  currentPass, debugConnectivityMatrix)

    if (debugConnectivityMatrix)
        checkConnectivityMatrix(obj.connectivityMatrix, targetSourceRFsNum, currentPass);
    end

    % Find all destinationRFs that have targetSourceRFsNum
    ss = squeeze(full(sum(obj.connectivityMatrix,1)));

    theTargetedDestinationRFindices = find(ss == targetSourceRFsNum);
    destinationRFsNumWithTargetSourceRFsNum = numel(theTargetedDestinationRFindices);
    if (destinationRFsNumWithTargetSourceRFsNum == 0)
        optimizationTransfersNum = 0;
        theOptimizedTargetDestinationRFindices = [];
        return;
    end

    % Sort theTargetedDestinationRFindices (based on their eccentricity & optimization center)
    idx = obj.sortDestinationRFsBasedOnOptimizationCenter(theTargetedDestinationRFindices);
    theTargetedDestinationRFindices = theTargetedDestinationRFindices(idx);

    % Keep track of which destination RFs were optimized
    theOptimizedTargetDestinationRFindices = [];
    optimizationTransfersNum = 0;

    % Go through each of theTargetedDestinationRFindices
    for iDestinationRFindex = 1:numel(theTargetedDestinationRFindices)
        theTargetDestinationRF = theTargetedDestinationRFindices(iDestinationRFindex);

        % Reset the involved destination and source RFs
        allInvolvedDestinationRFindices = [];
        allInvolvedSourceRFindices = [];

        % Add theTargetDestinationRF and its sourceRFs
        allInvolvedDestinationRFindices(1) = theTargetDestinationRF;
        allInvolvedSourceRFindices(1:targetSourceRFsNum) = find(squeeze(obj.connectivityMatrix(:,theTargetDestinationRF))>0.001);

        % Find the indices of the neigboring destinationRFs. These come sorted in increasing distance from theTargetDestinationRF
        [theNearbyDestinationRFindicesToTheTargetDestinationRF, theDistancesToTheTargetDestinationRF] = ...
            obj.indicesOfNeighboringDestinationRFs(...
                theTargetDestinationRF, ...
                'maxNeighborsNum', 6);

        % Accumulate valid nearbyDestinationRFs and their sourceRFs
        for iNearbyDestinationRF = 1:numel(theNearbyDestinationRFindicesToTheTargetDestinationRF)

            % Find the input cone indices for this nearbyDestinationRF
            theNearbyDestinationRFindex = theNearbyDestinationRFindicesToTheTargetDestinationRF(iNearbyDestinationRF);
            theNearbyDestinationRFinputSourceRFindices = find(squeeze(obj.connectivityMatrix(:,theNearbyDestinationRFindex))>0.001);

            if (minNumberOfSourceInputsToConsideredNearbyDestinationRFs > 0)
                % The nearby destinationRF must have the desired min number of inputs (at least targetSourceRFsNum+2)
                % so we can get one sourceRF from it
                if (numel(theNearbyDestinationRFinputSourceRFindices) < minNumberOfSourceInputsToConsideredNearbyDestinationRFs)
                    % Not enough inputs, dont include
                    continue;
                end
            end
            
            % Enough inputs or minNumberOfSourceInputsToConsideredNearbyDestinationRFs = 0, -> include
            allInvolvedDestinationRFindices = cat(1, allInvolvedDestinationRFindices(:), theNearbyDestinationRFindex);
            allInvolvedSourceRFindices = cat(1, allInvolvedSourceRFindices(:), theNearbyDestinationRFinputSourceRFindices(:));
        end % iNearbyDestinationRF

        % If we have more than 1 involved RGC try to optimize the involved connectivity matrix
        if (numel(allInvolvedDestinationRFindices) == 1)
            continue;
        end

        if (debugConnectivityMatrix)
            checkConnectivityMatrix(obj.connectivityMatrix(allInvolvedSourceRFindices, ...
                allInvolvedDestinationRFindices), targetSourceRFsNum, -10);
        end

        % Radius for searching for 
        searchRadius = 1.25*median(obj.destinationRFspacingsFromCentroids(theNearbyDestinationRFindicesToTheTargetDestinationRF));

        % Optimize the rgcRFcenterConeConnectivityMatrix for allInvolvedDestinationRFindices and allInvolvedSourceRFindices
        if (minNumberOfSourceInputsToConsideredNearbyDestinationRFs > 0)
            [updatedConnectivityMatrix, theDonorDestinationRFindex] = optimizeLocalConnectivityMatrixForSourceRFtransfer(...
                obj.connectivityMatrix(allInvolvedSourceRFindices, allInvolvedDestinationRFindices), ...
                obj.wiringParams.spatialChromaticUniformityTradeoff, ...
                obj.sourceLattice.metaData.coneTypes(allInvolvedSourceRFindices), ...
                obj.sourceLattice.RFpositionsMicrons(allInvolvedSourceRFindices,:), ...
                targetSourceRFsNum, searchRadius, debugConnectivityMatrix, ...
                obj.destinationRFspacingsFromCentroids(theTargetDestinationRF));
        else
            [updatedConnectivityMatrix, theDonorDestinationRFindex, costBeforeSwap, costAfterSwap] = optimizeLocalConnectivityMatrixForSourceRFswapping(...
                obj.connectivityMatrix(allInvolvedSourceRFindices, allInvolvedDestinationRFindices), ...
                obj.wiringParams.spatialChromaticUniformityTradeoff, ...
                obj.sourceLattice.metaData.coneTypes(allInvolvedSourceRFindices), ...
                obj.sourceLattice.RFpositionsMicrons(allInvolvedSourceRFindices,:), ...
                targetSourceRFsNum, searchRadius, debugConnectivityMatrix, ...
                obj.destinationRFspacingsFromCentroids(theTargetDestinationRF));
        end

        if (~isempty(updatedConnectivityMatrix))
            if (minNumberOfSourceInputsToConsideredNearbyDestinationRFs > 0) || ...
               ((minNumberOfSourceInputsToConsideredNearbyDestinationRFs == 0) && (costAfterSwap < costBeforeSwap))
            
                optimizationTransfersNum = optimizationTransfersNum + 1;
                obj.connectivityMatrix(allInvolvedSourceRFindices, allInvolvedDestinationRFindices) = updatedConnectivityMatrix;

                if (debugConnectivityMatrix)
                    checkConnectivityMatrix(obj.connectivityMatrix, targetSourceRFsNum, currentPass);
                end

                % Update the centroids of the involved destination RF indices
                theDestinationRFindicesInNeedOfCentroidUpdate = allInvolvedDestinationRFindices([1 theDonorDestinationRFindex]);
                obj.updateDestinationCentroidsFromInputs(theDestinationRFindicesInNeedOfCentroidUpdate);

                % Keep track of optimized destination RFindices
                theOptimizedTargetDestinationRFindices(numel(theOptimizedTargetDestinationRFindices)+1) = allInvolvedDestinationRFindices(1);
            end
        end
    end % for iDestinationRFindex

    % Update the destinationRF spacings based on the updated connectivity
    obj.updateDestinationRFspacingsBasedOnCentroids();
end

function [theUpdatedLocalConnectivityMatrix, theParentDestinationRFindexOfTheOptimalSourceRF] = optimizeLocalConnectivityMatrixForSourceRFtransfer(theLocalConnectivityMatrix, ...
    spatialChromaticUniformityTradeoff, sourceRFtypes, ...
    sourceRFpositions, targetSourceRFsNum, ...
    searchRadius, debugConnectivityMatrix, ...
    theTargetDestinationRFspacing)

    theUpdatedLocalConnectivityMatrix = [];
    theParentDestinationRFindexOfTheOptimalSourceRF = [];

    % We will compute costs for adding each and every of the surrounding sourceRFs
    theSurroundingSourceRFindices = (targetSourceRFsNum+1):size(sourceRFpositions,1);
    theSurroundingSourceRFpositions = sourceRFpositions(theSurroundingSourceRFindices,:);
    
    % The targetDestinationRF is always at index1
    theTargetDestinationRFindex = 1;

    % The source RFs for theTargetDestinationRFindex are always the first targetSourceRFsNum
    theTargetDestinationRFsourceRFindices = 1:targetSourceRFsNum;
    theTargetDestinationRFcentroid = mean(sourceRFpositions(theTargetDestinationRFsourceRFindices,:),1);

    % Examine each and every of theSurroundingSourceRFindices
    theConsideredParentDestinationRFindices = []; 
    theConsideredSourceRFindices = []; 
    theConsideredCosts = [];

    for idx = 1:numel(theSurroundingSourceRFindices)
        % Only examine surrounding sourceRFs whose distance to theTargetDestinationRFcentroid is < searchRadius
        surroundingSourceRFdistancesToTargetCentroid = ...
            sqrt(sum((bsxfun(@minus, theSurroundingSourceRFpositions(idx,:), theTargetDestinationRFcentroid)).^2,2));

        if (surroundingSourceRFdistancesToTargetCentroid > searchRadius)
            % Too far from targetRF centroid. Exclude this one.
            continue;
        end

        % Find which destinationRF this candidate sourceRF is connected to
        iNearbyCandidateSourceRFindex = theSurroundingSourceRFindices(idx);
        theParentDestinationRFindex = find(full(squeeze(theLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, :))) > 0.001);

        if (isempty(theParentDestinationRFindex))
            % Not connected to a destinationRF. Probaby a surround cone. Exclude this one.
            continue;
        end

        % Check that theParentDestinationRFindex has at least 2+ source RFs than theTargetDestinationRFindex
        theSourceRFsNum = find(full(squeeze(theLocalConnectivityMatrix(:,theTargetDestinationRFindex))) > 0.001);
        if (theSourceRFsNum ~= targetSourceRFsNum)
            error('Error in logic');
        end
        theParentSourceRFsNum = find(full(squeeze(theLocalConnectivityMatrix(:,theParentDestinationRFindex))) > 0.001);
        if (~(theParentSourceRFsNum >= targetSourceRFsNum+2))
            fprintf('Yeap, this would have been an invalid transfer\n');
            % theParentDestinationRFindex does not have at least 2+ source RFs than theTargetDestinationRFindex
            continue;
        end

        if (theParentDestinationRFindex == theTargetDestinationRFindex)
            % Big problem with logic
            error('oh no')
        end

        % Assemble the candidate connectivity matrix to compute the cost of transfering
        % iNearbyCandidateSourceRFindex from theParentDestinationRFindex to theTargetDestinationRFindex
        theCandidateLocalConnectivityMatrix = theLocalConnectivityMatrix;
        % Disconnect iNearbyCandidateSourceRFindex from its original parent
        theCandidateLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, theParentDestinationRFindex) = 0;
        % Connect iNearbyCandidateSourceRFindex to theTargetDestinationRFindex 
        theCandidateLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, theTargetDestinationRFindex) = 1;
        
        % Compute costs for this candidate tranfer
        [theConsideredCosts(numel(theConsideredCosts )+1), theSpatialCompactnessCost, theSourceTypePurityCost, theCentroidOverlapCost, theVarianceCost] = ...
            costForConnectivityMatrix(theCandidateLocalConnectivityMatrix, spatialChromaticUniformityTradeoff, ...
                            sourceRFpositions, sourceRFtypes, theParentDestinationRFindex, ...
                            theTargetDestinationRFspacing, theTargetDestinationRFspacing);

        theConsideredSourceRFindices(numel(theConsideredSourceRFindices)+1) = iNearbyCandidateSourceRFindex;
        theConsideredParentDestinationRFindices(numel(theConsideredParentDestinationRFindices)+1) = theParentDestinationRFindex;
    end % for idx

    if (isempty(theConsideredCosts))
        return;
    end

    % Determine theOptimalSourceRFindex
    [minCost,idx] = min(theConsideredCosts);
    if (isinf(minCost))
        return;
    end

    theOptimalSourceRFindex = theConsideredSourceRFindices(idx);
    theParentDestinationRFindexOfTheOptimalSourceRF = theConsideredParentDestinationRFindices(idx);

    % Assemble the updated connectivity matrix 
    theUpdatedLocalConnectivityMatrix = theLocalConnectivityMatrix;
    % Disconnect theOptimalSourceRFindex from its original parent
    theUpdatedLocalConnectivityMatrix(theOptimalSourceRFindex, theParentDestinationRFindexOfTheOptimalSourceRF) = 0;
    % Connect theOptimalSourceRFindex to theTargetDestinationRFindex 
    theUpdatedLocalConnectivityMatrix(theOptimalSourceRFindex, theTargetDestinationRFindex) = 1;

    if (debugConnectivityMatrix)
        checkConnectivityMatrix(theLocalConnectivityMatrix, targetSourceRFsNum, -1);
        checkConnectivityMatrix(theUpdatedLocalConnectivityMatrix, targetSourceRFsNum, -2);
    end
end


function [theUpdatedLocalConnectivityMatrix, theParentDestinationRFindexOfTheOptimalSourceRF, costBeforeSwap, costAfterSwap] = ...
        optimizeLocalConnectivityMatrixForSourceRFswapping(...
            theLocalConnectivityMatrix, spatialChromaticUniformityTradeoff, sourceRFtypes, sourceRFpositions, targetSourceRFsNum, ...
            searchRadius, debugConnectivityMatrix, theTargetDestinationRFspacing)

    theUpdatedLocalConnectivityMatrix = [];
    theParentDestinationRFindexOfTheOptimalSourceRF = [];
     
    theSurroundingSourceRFindices = (targetSourceRFsNum+1):size(sourceRFpositions,1);
    theSurroundingSourceRFpositions = sourceRFpositions(theSurroundingSourceRFindices,:);

    % The targetDestinationRF is always at index1
    theTargetDestinationRFindex = 1;
    
    theTargetDestinationRFsourceRFindices = 1:targetSourceRFsNum;
    theTargetDestinationRFsourceRFpositions = sourceRFpositions(theTargetDestinationRFsourceRFindices,:);

    % Examine each and every combination of [targetSourceRFs x surrounding sourceRFs]
    theConsideredParentDestinationRFindices = [];  
    theConsideredSourceRFindices = []; 
    theConsideredTargetDestinationRFsourceRFindices = []; 
    theConsideredCosts = [];
    theBeforeSwapCosts = [];

    for idxSurround = 1:numel(theSurroundingSourceRFindices)
        % Find which destinationRF this candidate sourceRF is connected to
        iNearbyCandidateSourceRFindex = theSurroundingSourceRFindices(idxSurround);
        theParentDestinationRFindex = find(full(squeeze(theLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, :))) > 0.001);

        if (isempty(theParentDestinationRFindex))
            % Not connected to a destinationRF. Probaby a surround cone. Exclude this one.
            continue;
        end

        if (theParentDestinationRFindex == theTargetDestinationRFindex)
            % Big problem with logic
            error('oh no')
        end

        for idxTarget = 1:numel(theTargetDestinationRFsourceRFindices)

            % Compute the before swap cost
            theBeforeSwapCosts(numel(theBeforeSwapCosts)+1) = costForConnectivityMatrix(...
                    theLocalConnectivityMatrix, spatialChromaticUniformityTradeoff, ...
                    sourceRFpositions, sourceRFtypes, theParentDestinationRFindex, ...
                    theTargetDestinationRFspacing, theTargetDestinationRFspacing);

            iTargetCandidateSourceRFindex = theTargetDestinationRFsourceRFindices(idxTarget);
            % Assemble the candidate connectivity matrix to compute the cost of swapping
            % iNearbyCandidateSourceRFindex from theParentDestinationRFindex 
            % with
            % iTargetCandidateSourceRFindex from theTargetDestinationRFindex 

            theCandidateLocalConnectivityMatrix = theLocalConnectivityMatrix;
            % Disconnect iNearbyCandidateSourceRFindex from its original parent
            theCandidateLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, theParentDestinationRFindex) = 0;
            % Connect iNearbyCandidateSourceRFindex to theTargetDestinationRFindex 
            theCandidateLocalConnectivityMatrix(iNearbyCandidateSourceRFindex, theTargetDestinationRFindex) = 1;

            % Disconnect iTargetCandidateSourceRFindex from theTargetDestinationRFindex 
            theCandidateLocalConnectivityMatrix(iTargetCandidateSourceRFindex, theTargetDestinationRFindex) = 0;
            % And connect iTargetCandidateSourceRFindex to theParentDestinationRFindex
            theCandidateLocalConnectivityMatrix(iTargetCandidateSourceRFindex, theParentDestinationRFindex) = 1;

            % Compute costs for this candidate swap
            [theConsideredCosts(numel(theConsideredCosts )+1), theSpatialCompactnessCost, theSourceTypePurityCost] = ...
                costForConnectivityMatrix(theCandidateLocalConnectivityMatrix, spatialChromaticUniformityTradeoff, ...
                            sourceRFpositions, sourceRFtypes, theParentDestinationRFindex, ...
                            theTargetDestinationRFspacing, theTargetDestinationRFspacing);

            theConsideredTargetDestinationRFsourceRFindices(numel(theConsideredTargetDestinationRFsourceRFindices)+1) = ...
                iTargetCandidateSourceRFindex;
            theConsideredSourceRFindices(numel(theConsideredSourceRFindices)+1) = ...
                iNearbyCandidateSourceRFindex;
            theConsideredParentDestinationRFindices(numel(theConsideredParentDestinationRFindices)+1) = ...
                theParentDestinationRFindex;
        end % idxTarget
    end % idxSurround

    if (isempty(theConsideredCosts))
        return;
    end

    % Determine the optimal combination
    [minCost,idx] = min(theConsideredCosts);
    if (isinf(minCost))
        return;
    end

    theOptimalTargetDestinationRFsourceRFindex = theConsideredTargetDestinationRFsourceRFindices(idx);
    theOptimalSourceRFindex = theConsideredSourceRFindices(idx);
    theParentDestinationRFindexOfTheOptimalSourceRF = theConsideredParentDestinationRFindices(idx);

    costBeforeSwap = theBeforeSwapCosts(idx);
    costAfterSwap = theConsideredCosts(idx);

    % Assemble the updated connectivity matrix 
    theUpdatedLocalConnectivityMatrix = theLocalConnectivityMatrix;

    % Disconnect theOptimalSourceRFindex from its original parent
    theUpdatedLocalConnectivityMatrix(theOptimalSourceRFindex, theParentDestinationRFindexOfTheOptimalSourceRF) = 0;
    % Connect theOptimalSourceRFindex to theTargetDestinationRFindex 
    theUpdatedLocalConnectivityMatrix(theOptimalSourceRFindex, theTargetDestinationRFindex) = 1;

    % Disconnect theOptimalTargetDestinationRFsourceRFindex from theTargetDestinationRFindex 
    theUpdatedLocalConnectivityMatrix(theOptimalTargetDestinationRFsourceRFindex, theTargetDestinationRFindex ) = 0;
    % And connect theOptimalTargetDestinationRFsourceRFindex to theParentDestinationRFindex
    theUpdatedLocalConnectivityMatrix(theOptimalTargetDestinationRFsourceRFindex, theParentDestinationRFindexOfTheOptimalSourceRF) = 1;

    if (debugConnectivityMatrix)
        checkConnectivityMatrix(theLocalConnectivityMatrix, targetSourceRFsNum, -1);
        checkConnectivityMatrix(theUpdatedLocalConnectivityMatrix, targetSourceRFsNum, -2);
    end
end


function [theTotalCost, theSpatialCompactnessCost, theSourceTypeUniformityCost, theCentroidOverlapCost, theVarianceCost] = ...
        costForConnectivityMatrix(theLocalConnectivityMatrix, spatialChromaticUniformityTradeoff, ...
                        sourceRFpositions, sourceRFtypes, theDonorDestinationRFindex, ...
                        theTargetDestinationRFspacing, theDonorDestinationRFspacing)

    % Determine positions and types of sourceRFs for theTargetDestinationRFindex 
    theTargetDestinationRFindex = 1;
    idx = find(full(squeeze(theLocalConnectivityMatrix(:, theTargetDestinationRFindex)))>0.001);
    theTargetRFsourceRFtypes = sourceRFtypes(idx);
    theTargetRFSourceRFpositions = sourceRFpositions(idx,:);

    % Determine positions and types of sourceRFs for theDonorDestinationRFindex
    idx = find(full(squeeze(theLocalConnectivityMatrix(:, theDonorDestinationRFindex)))>0.001);
    theDonorRFsourceRFtypes = sourceRFtypes(idx);
    theDonorRFSourceRFpositions = sourceRFpositions(idx,:);

    % Compute the theSpatialCompactnessCost for this pair of destinationRFs
    [theSpatialCompactnessCost, ~, theCentroidOverlapCost, theVarianceCost] = coneToMidgetRGCConnector.spatialCompactnessCost(...
        theTargetRFSourceRFpositions, theDonorRFSourceRFpositions, [], [], ...
        theTargetDestinationRFspacing, theDonorDestinationRFspacing);

    % Compute theSourceTypeUniformityCost for this pair of destinationRFs
    theTargetSourceTypeUniformityCost = coneToMidgetRGCConnector.spectralUniformityCost(theTargetRFsourceRFtypes, []);
    theDonorSourceTypeUniformityCost = coneToMidgetRGCConnector.spectralUniformityCost(theDonorRFsourceRFtypes, []);
    theSourceTypeUniformityCost = 0.5*(theTargetSourceTypeUniformityCost+theDonorSourceTypeUniformityCost);

    % The total cost: weighted sum of theSpatialCompactnessCost and theSourceTypePurityCost
    theTotalCost = ...
        spatialChromaticUniformityTradeoff     * theSpatialCompactnessCost  + ...
        (1-spatialChromaticUniformityTradeoff) * theSourceTypeUniformityCost;
end