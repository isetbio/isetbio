function transferConesBetweenNearbyRGCsWithUnbalancedInputNumerosities(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('optimizationCenter', 'patchCenter', @(x)(ismember(x, {'patchCenter', 'visualFieldCenter'})));

    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;
    optimizationCenter = p.Results.optimizationCenter;

    minConeInputsNum = 3;
    maxConeInputsNum = 6;

    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step3', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end


    currentPass = 0;
    tranfersInCurrentPass = 1;
    convergenceAchieved = false;

    % Compute the current costs for all RGCs to maintain their current cone inputs 
    [totalCosts, spatialCosts, chromaticCosts] = obj.computeInputMaintenanceCostAcrossEntireMosaic();

    % Net cost across all RGCs (Initial)
    netTotalCostInitial = mean(totalCosts);
    netSpatialCostInitial = mean(spatialCosts(spatialCosts>=0));
    netChromaticCostInitial = mean(chromaticCosts(chromaticCosts>=0));

    while (currentPass < obj.wiringParams.maxPassesNum) && (tranfersInCurrentPass>0) && (~convergenceAchieved)

        % Update pass number
        currentPass = currentPass + 1;
        tranfersInCurrentPass = 0;

        % We only consider RGCs that have between [minConeInputsNum maxConeInputsNum] input cones
        ss = squeeze(sum(obj.coneConnectivityMatrix,1));
        targetedRGCindices = find((ss >= minConeInputsNum) & (ss <= maxConeInputsNum));
        if (isempty(targetedRGCindices))
            continue;
        end

        % Retrieve the centroids of the targered RGCs
        targetedRGCCentroids = obj.RGCRFcentroidsFromInputs(targetedRGCindices,:);
    
        % Sort targeted RGCs according to their ecc
        switch (optimizationCenter)
            case 'visualFieldCenter'
                ecc = sum(targetedRGCCentroids.^2,2);
            case 'patchCenter'
                diff = bsxfun(@minus, targetedRGCCentroids, obj.inputConeMosaic.eccentricityMicrons);
                ecc = sum(diff.^2,2);
        end
        [~, idx] = sort(ecc, 'ascend');
        sortedRGCindices = targetedRGCindices(idx);
        
    
        for iRGC = 1:numel(sortedRGCindices)
            theSourceRGCindex = sortedRGCindices(iRGC);
            
            % Find up to N neigboring RGCs that are no farther than k x RGC separation
            nearbyRGCindices = obj.neihboringRGCindices(theSourceRGCindex);
            if (isempty(nearbyRGCindices))
                continue;
            end

            % Find the indices and weights of cone inputs to this RGC
            theSourceRGCinputConeIndices = find(squeeze(obj.coneConnectivityMatrix(:,theSourceRGCindex))>0);
            theSourceRGCinputConesNum = numel(theSourceRGCinputConeIndices);
            theSourceRGCinputConeWeights = full(obj.coneConnectivityMatrix(theSourceRGCinputConeIndices,theSourceRGCindex));
    
            

            % Find the difference in # of input cones between source RGC and neigboring RGCs
            inputConeIndicesAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
            inputConeWeightsAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
            differenceInConeInputs = zeros(1,numel(nearbyRGCindices));
            for iNearbyRGC = 1:numel(nearbyRGCindices)
                theNearbyRGCindex = nearbyRGCindices(iNearbyRGC);
                inputConeIndicesAllNeighboringRGCs{iNearbyRGC} = find(obj.coneConnectivityMatrix(:,theNearbyRGCindex)>0);
                inputConeWeightsAllNeighboringRGCs{iNearbyRGC} = full(obj.coneConnectivityMatrix(inputConeIndicesAllNeighboringRGCs{iNearbyRGC},theNearbyRGCindex));
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
                obj.optimizeTransferOfConeInputs(...
                    theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
                    theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivityState(1002);
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end %if (maxDiff >= 2)
        end % iRGC

        % Update the local RGCRFspacings based on updated connectivity
        obj.updateLocalRGCRFspacingsBasedOnCurrentCentroids();

        % Compute the current costs for all RGCs to maintain their current cone inputs 
        [totalCosts, spatialCosts, chromaticCosts] = obj.computeInputMaintenanceCostAcrossEntireMosaic();

        netTransfers(currentPass) = tranfersInCurrentPass;

        % Net cost across all RGCs
        netTotalCost(currentPass) = mean(totalCosts);
        netSpatialCost(currentPass) = mean(spatialCosts(spatialCosts>=0));
        netChromaticCost(currentPass) = mean(chromaticCosts(chromaticCosts>=0));
        
        % Check delta in netTotalCost for convergence
        if (currentPass == 1)
            lastTotalCost = netTotalCostInitial;
        else
            lastTotalCost = netTotalCost(currentPass-1);
        end

        netTotalCostSequence = [netTotalCostInitial netTotalCost];


        
        % Visualize convergence
        RGCconnector.visualizeConvergence(currentPass, netTotalCostInitial, netTotalCost, ...
                                          netSpatialCostInitial, netSpatialCost, ...
                                          netChromaticCostInitial, netChromaticCost, ...
                                          netTransfers, obj.wiringParams.maxPassesNum);

        
        % Determine whether convergence was achieved
        netTotalCostSequence = [netTotalCostInitial netTotalCost];
        convergenceAchieved = RGCconnector.convergenceAchieved(netTotalCostSequence);

    end % while loop

    if (generateProgressVideo)
        videoOBJ.close();
    end

end