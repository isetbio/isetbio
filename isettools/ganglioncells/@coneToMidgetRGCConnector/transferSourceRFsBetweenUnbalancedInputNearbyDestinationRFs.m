function transferSourceRFsBetweenUnbalancedInputNearbyDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;

    % We only transfer cones from RGC-A to a nearby RGC-B, if RGC-A has between 3 and 6 cones
    % If it has less than 3, i.e. 2, then this will be dealt by the
    % swapSourceRFs method which is called later on
    minConeInputsNum = 3;
    maxConeInputsNum = 6;

    ss = squeeze(sum(obj.connectivityMatrix,1));
    targetedRGCindices = find((ss >= minConeInputsNum) & (ss <= maxConeInputsNum))
    if (isempty(targetedRGCindices))
        fprintf('No RGCs with 3-6 input cones. Transfer cone inputs phase will not execute. \n.')
        return;
    end

    if (strcmp(obj.wiringParams.optimizationCenter, 'patchCenter'))
        coneMosaicCenter = mean(obj.sourceLattice.RFpositionsMicrons,1);
    end

    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step3', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    % Compute the current costs for all RGCs to maintain their current cone inputs 
    theCostComponents = obj.totalInputMaintenanceCost();
    totalCosts = theCostComponents(:,1);
    spatialCosts = theCostComponents(:,2);
    chromaticCosts = theCostComponents(:,3);

    % Net cost across all RGCs (Initial)
    netTotalCostInitial = mean(totalCosts);
    netSpatialCostInitial = mean(spatialCosts(spatialCosts>=0));
    netChromaticCostInitial = mean(chromaticCosts(chromaticCosts>=0));


    currentPass = 0;  tranfersInCurrentPass = 1; convergenceAchieved = false;
    while (currentPass < obj.wiringParams.maxPassesNum) && ...
          (tranfersInCurrentPass>0) && ...
          (~convergenceAchieved)

        % Update pass number
        currentPass = currentPass + 1;
        tranfersInCurrentPass = 0;


        % We only consider RGCs that have between [minConeInputsNum maxConeInputsNum] input cones
        ss = squeeze(sum(obj.connectivityMatrix,1));
        targetedRGCindices = find((ss >= minConeInputsNum) & (ss <= maxConeInputsNum));
        if (isempty(targetedRGCindices))
            continue;
        end

        % Retrieve the centroids of the targered RGCs
        targetedRGCCentroids = obj.destinationRFcentroidsFromInputs(targetedRGCindices,:);
    
        switch (obj.wiringParams.optimizationCenter)
            case 'visualFieldCenter'
                ecc = sum(targetedRGCCentroids.^2,2);
            case 'patchCenter'
                diff = bsxfun(@minus, targetedRGCCentroids, coneMosaicCenter);
                ecc = sum(diff.^2,2);
        end % switch

        [~, idx] = sort(ecc, 'ascend');
        sortedRGCindices = targetedRGCindices(idx);

        for iRGC = 1:numel(sortedRGCindices)
            theSourceRGCindex = sortedRGCindices(iRGC);
            
            % Find the indicies of the neigboring RGCs
            nearbyRGCindices = obj.indicesOfNeighboringDestinationRFs(theSourceRGCindex);
            if (isempty(nearbyRGCindices))
                continue;
            end

            % Find the indices and weights of cone inputs to this RGC
            theSourceRGCinputConeIndices = find(squeeze(obj.connectivityMatrix(:,theSourceRGCindex))>0);
            theSourceRGCinputConesNum = numel(theSourceRGCinputConeIndices);
            theSourceRGCinputConeWeights = full(obj.connectivityMatrix(theSourceRGCinputConeIndices,theSourceRGCindex));
    
            % Find the difference in # of input cones between source RGC and neigboring RGCs
            inputConeIndicesAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
            inputConeWeightsAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
            differenceInConeInputs = zeros(1,numel(nearbyRGCindices));
            for iNearbyRGC = 1:numel(nearbyRGCindices)
                theNearbyRGCindex = nearbyRGCindices(iNearbyRGC);
                inputConeIndicesAllNeighboringRGCs{iNearbyRGC} = find(obj.connectivityMatrix(:,theNearbyRGCindex)>0);
                inputConeWeightsAllNeighboringRGCs{iNearbyRGC} = full(obj.connectivityMatrix(inputConeIndicesAllNeighboringRGCs{iNearbyRGC},theNearbyRGCindex));
                differenceInConeInputs(iNearbyRGC) = theSourceRGCinputConesNum - numel(inputConeIndicesAllNeighboringRGCs{iNearbyRGC});
            end % iNearbyRGC

            % Find nearby RGCs whose # of input cones is at least 2 less
            maxDiff = max(differenceInConeInputs);
            if (maxDiff >= 2)
                tranfersInCurrentPass = tranfersInCurrentPass + 1;
                idx = find(differenceInConeInputs == maxDiff);
                allNeighboringRGCsInputConeIndices = inputConeIndicesAllNeighboringRGCs(idx);
                allNeighboringRGCsInputConeWeights = inputConeWeightsAllNeighboringRGCs(idx);
                theNeighboringRGCindices = nearbyRGCindices(idx);
                obj.optimizeTransferOfInputRFs(...
                    theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
                    theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivity(1002);
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end %if (maxDiff >= 2)
        end % for iRGC

        % Update the RF spacings based on the updated connectivity
        obj.updateDestinationRFspacingsBasedOnCentroids();

        
        % Compute the current costs for all RGCs to maintain their current cone inputs 
        theCostComponents = obj.totalInputMaintenanceCost();
        totalCosts = theCostComponents(:,1);
        spatialCosts = theCostComponents(:,2);
        chromaticCosts = theCostComponents(:,3);

        % Net cost across all RGCs (Initial)
        netTotalCost(currentPass) = mean(totalCosts);
        netSpatialCost(currentPass) = mean(spatialCosts(spatialCosts>=0));
        netChromaticCost(currentPass) = mean(chromaticCosts(chromaticCosts>=0));

        netTransfers(currentPass) = tranfersInCurrentPass;

        % Check delta in netTotalCost for convergence
        if (currentPass == 1)
            lastTotalCost = netTotalCostInitial;
        else
            lastTotalCost = netTotalCost(currentPass-1);
        end

        % Accumulate cost sequence
        netTotalCostSequence = [netTotalCostInitial netTotalCost];
        netSpatialCostSequence = [netSpatialCostInitial netSpatialCost];
        netChromaticCostSequence = [netChromaticCostInitial netChromaticCost];
        
        % Visualize convergence
        costsMatrix(:,1) = netTotalCostSequence(:);
        costsMatrix(:,2) = netSpatialCostSequence(:);
        costsMatrix(:,3) = netChromaticCostSequence(:);
        costsNames = {'total cost', 'spatial variance cost', 'chromatic variance cost'};
        MosaicConnector.visualizeConvergenceSequence(currentPass, ...
            costsMatrix, costsNames, ...
            netTransfers, obj.wiringParams.maxPassesNum);


        % Determine whether convergence was achieved
        convergenceAchieved = MosaicConnector.convergenceAchieved(netTotalCostSequence);

    end % while loop

    if (generateProgressVideo)
        videoOBJ.close();
    end

end
