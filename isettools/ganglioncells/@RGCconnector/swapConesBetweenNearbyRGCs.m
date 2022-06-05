function swapConesBetweenNearbyRGCs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('optimizationCenter', 'patchCenter', @(x)(ismember(x, {'patchCenter', 'visualFieldCenter'})));

    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;
    optimizationCenter = p.Results.optimizationCenter;

    % Video setup
    
    if (generateProgressVideo)
        videoOBJ = VideoWriter(sprintf('Step5_w%2.3f', obj.wiringParams.chromaticSpatialVarianceTradeoff), 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    else
        videoOBJ = [];
    end


    currentPass = 0;
    swapsInCurrentPass = 1;
    convergenceAchieved = false;

    [totalCosts, spatialCosts, chromaticCosts] = obj.computeInputMaintenanceCostAcrossEntireMosaic();

    % Net cost across all RGCs (Initial)
    netTotalCostInitial = mean(totalCosts);
    netSpatialCostInitial = mean(spatialCosts(spatialCosts>=0));
    netChromaticCostInitial = mean(chromaticCosts(chromaticCosts>=0));
           

    while (currentPass < obj.wiringParams.maxPassesNum) && (swapsInCurrentPass>0) && (~convergenceAchieved)
        % Update pass number
        currentPass = currentPass + 1;
        swapsInCurrentPass = 0;

        % Retrieve the current centroids of all RGCs
        targetedRGCCentroids = obj.RGCRFcentroidsFromInputs;
    
        % Sort RGCs according to their ecc
        switch (optimizationCenter)
            case 'visualFieldCenter'
                ecc = sum(targetedRGCCentroids.^2,2);
            case 'patchCenter'
                diff = bsxfun(@minus, targetedRGCCentroids, obj.inputConeMosaic.eccentricityMicrons);
                ecc = sum(diff.^2,2);
        end
        [~,sortedRGCindices] = sort(ecc, 'ascend');


        for iRGC = 1:numel(sortedRGCindices)
            fprintf('[cone swapping phase (pass #%d)]:: evaluating benefits for RGC %d/%d\n', currentPass, iRGC, numel(sortedRGCindices));

            % Find up to N neigboring RGCs that are no farther than k x RGC separation
            theSourceRGCindex = sortedRGCindices(iRGC);
            theNeighboringRGCindices = obj.neihboringRGCindices(theSourceRGCindex);
            if (isempty(theNeighboringRGCindices))
                continue;
            end

            % Find the indices and weights of cone inputs to this RGC
            theSourceRGCinputConeIndices = find(squeeze(obj.coneConnectivityMatrix(:,theSourceRGCindex))>0);
            theSourceRGCinputConeWeights = full(obj.coneConnectivityMatrix(theSourceRGCinputConeIndices,theSourceRGCindex));
    
            if (numel(theSourceRGCinputConeIndices)==1)
                % Dont do anything if the source has a single cone input
                continue;
            end

            % Find the indices and weights of cone inputs to all neigboring RGCs
            allNeighboringRGCsInputConeIndices = cell(1,numel(theNeighboringRGCindices));
            allNeighboringRGCsInputConeWeights= cell(1,numel(theNeighboringRGCindices));

            % Find the indices of input cones to all neigboring RGCs. 
            % Also compute the mean # of cone inputs near this RGC
            meanConeInputsNum = numel(theSourceRGCinputConeIndices);
            for iNearbyRGC = 1:numel(theNeighboringRGCindices)
                theNearbyRGCindex = theNeighboringRGCindices(iNearbyRGC);
                allNeighboringRGCsInputConeIndices{iNearbyRGC} = find(obj.coneConnectivityMatrix(:,theNearbyRGCindex)>0);
                allNeighboringRGCsInputConeWeights{iNearbyRGC} = full(obj.coneConnectivityMatrix(allNeighboringRGCsInputConeIndices{iNearbyRGC},theNearbyRGCindex));
                meanConeInputsNum = meanConeInputsNum + numel(allNeighboringRGCsInputConeIndices{iNearbyRGC});
            end % iNearbyRGC
            meanConeInputsNum = meanConeInputsNum / (1+numel(theNeighboringRGCindices));
            
            % Check whether the mean # of cone inputs < obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping)
            if (meanConeInputsNum > obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping)
                fprintf('No swapping for RGC %d of %d. Neighborhood RGCs have an average of %2.1f cone inputs. Max for swapping: %d.\n', ...
                    iRGC, numel(sortedRGCindices), meanConeInputsNum, obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping);
                continue;
            end

            % OK, mean # of cone inputs not too large, examine if there is a beneficial swap
            beneficialSwapWasFound = obj.optimizeSwappingOfConeInputs(...
                    theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
                    theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

            
            if (beneficialSwapWasFound)
                swapsInCurrentPass = swapsInCurrentPass + 1;
                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivityState(1005, ...
                        'titleString', sprintf('Cone swapping phase (PASS:%d, RGC:%d/%d)', currentPass, iRGC, numel(sortedRGCindices)));
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end

        end % iRGC

        % Update the local RGCRFspacings based on updated connectivity
        obj.updateLocalRGCRFspacingsBasedOnCurrentCentroids();


        % Compute the current costs for all RGCs to maintain their current cone inputs
        [totalCosts, spatialCosts, chromaticCosts] = obj.computeInputMaintenanceCostAcrossEntireMosaic();
        netSwaps(currentPass) = swapsInCurrentPass;

        % Net cost across all RGCs
        netTotalCost(currentPass) = mean(totalCosts);
        netSpatialCost(currentPass) = mean(spatialCosts(spatialCosts>=0));
        netChromaticCost(currentPass) = mean(chromaticCosts(chromaticCosts>=0));
        

        % Visualize convergence
        RGCconnector.visualizeConvergence(currentPass, netTotalCostInitial, netTotalCost, ...
                                          netSpatialCostInitial, netSpatialCost, ...
                                          netChromaticCostInitial, netChromaticCost, ...
                                          netSwaps, obj.wiringParams.maxPassesNum);

        % Determine whether convergence was achieved
        netTotalCostSequence = [netTotalCostInitial netTotalCost];
        convergenceAchieved = RGCconnector.convergenceAchieved(netTotalCostSequence);

        %close(hProgressBar);
    end % while pass loop

    if (generateProgressVideo)
        videoOBJ.close();
    end
end