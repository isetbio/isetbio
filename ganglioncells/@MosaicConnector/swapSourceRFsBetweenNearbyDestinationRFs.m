function swapSourceRFsBetweenNearbyDestinationRFs(obj, varargin)
     % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;


    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step5', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    % Compute the current costs for all destination RF to maintain their current inputs 
    % theCostComponents is an [nDestinationRFs x cost components]  matrix
    theCostComponentsMatrix = obj.totalInputMaintenanceCost();

    % The different cost component sequences
    s = theCostComponentsMatrix(1,:);
    idx = find(s > -99);
    netCostSequences = mean(theCostComponentsMatrix(idx,:),1);

    

    currentPass = 0;  swapsInCurrentPass = 1; convergenceAchieved = false;
    while (currentPass < obj.wiringParams.maxPassesNum) && (swapsInCurrentPass>0) && (~convergenceAchieved)

        % Save old connectivity matrix
        lastConnectivityMatrix = obj.connectivityMatrix;
        lastDestinationRFspacingsFromCentroids = obj.destinationRFspacingsFromCentroids;

        % Update pass number
        currentPass = currentPass + 1;
        swapsInCurrentPass = 0;

        allDestinationRFindices = 1:size(obj.destinationRFcentroidsFromInputs,1);
        sortedDestinationRFindices = obj.sortDestinationRFsBasedOnOptimizationCenter(allDestinationRFindices);

        % Feedback
        fprintf('\nEvaluating benefit of swapping inputs with nearby units for %d destination RFs (PASS %d/%d)\n',...
            numel(sortedDestinationRFindices),currentPass, obj.wiringParams.maxPassesNum);

        for iDestinationRF = 1:numel(sortedDestinationRFindices)
            % Feedback
            fprintf('.');
            if (mod((iDestinationRF-1),100) == 99)
                fprintf('  [%d/%d]\n', iDestinationRF , numel(sortedDestinationRFindices));
            end

            % Retrieve the destinationRF index 
            theSourceDestinationRFindex = sortedDestinationRFindices(iDestinationRF);
         
            % Find the indicies of the neigboring destinationRFs
            nearbyDestinationRFindices = obj.indicesOfNeighboringDestinationRFs(theSourceDestinationRFindex);
            if (isempty(nearbyDestinationRFindices))
                continue;
            end

            % Find the indices and weights of cone inputs to this destination RF
            theSourceDestinationRFinputIndices = find(squeeze(obj.connectivityMatrix(:,theSourceDestinationRFindex))>0);
            theSourceDestinationRFinputWeights = full(obj.connectivityMatrix(theSourceDestinationRFinputIndices,theSourceDestinationRFindex));
            % theSourceDestinationRFinputConesNum = numel(theSourceDestinationRFinputIndices);

            if (numel(theSourceDestinationRFinputIndices)==1)
                % Dont do anything if the destinationRF has a single input
                continue;
            end

            % Find the indices and weights of inputs to all neigboring RGCs
            inputIndicesAllNearbyDestinationRFs = cell(1,numel(nearbyDestinationRFindices));
            inputWeightsAllNearbyDestinationRFs = cell(1,numel(nearbyDestinationRFindices));

            % Find the indices of inputs to all nearby destination RFs
            % Also compute the mean # of inputs among the nearby destination RFs
            inputsNum = numel(nearbyDestinationRFindices);
            for iNearbyDestinationRF = 1:numel(nearbyDestinationRFindices)
                theNearbyDestinationRFindex = nearbyDestinationRFindices(iNearbyDestinationRF);
                inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF} = find(obj.connectivityMatrix(:,theNearbyDestinationRFindex)>0);
                inputWeightsAllNearbyDestinationRFs{iNearbyDestinationRF} = full(obj.connectivityMatrix(inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF},theNearbyDestinationRFindex));
                inputsNum(iNearbyDestinationRF) = numel(inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF});
            end % iNearbyDestinationRF
            meanInputsNum = floor(mean(inputsNum));

            % Check whether the mean # of inputs < obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping)
            if (meanInputsNum > obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping)
                fprintf('\nNo swapping for destination RF #%d Nearby destination RFs have an average of %2.1f inputs. Max for swapping: %d.\n', ...
                    iDestinationRF, meanInputsNum, obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping);
                continue;
            end
            
            % OK, mean # of cone inputs not too large, examine if there is a beneficial swap
            beneficialSwapWasFound = obj.optimizeSwappingOfInputRFs(...
                theSourceDestinationRFindex, theSourceDestinationRFinputIndices, theSourceDestinationRFinputWeights, ...
                nearbyDestinationRFindices, inputIndicesAllNearbyDestinationRFs, inputWeightsAllNearbyDestinationRFs);

            if (beneficialSwapWasFound)
                swapsInCurrentPass = swapsInCurrentPass + 1;
                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivity(1005, ...
                        'titleString', sprintf('Input swapping phase (PASS:%d, destinationRF:%d/%d)', ...
                        currentPass, iDestinationRF, numel(sortedDestinationRFindices)));
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end
        end % for iDestinationRF
        
        % Feedback
        fprintf('  [%d/%d]\n', iDestinationRF , numel(sortedDestinationRFindices));
        
        % Update the destinationRF spacings based on the updated connectivity
        obj.updateDestinationRFspacingsBasedOnCentroids();

        % Update swaps sequence
        netSwaps(currentPass) = swapsInCurrentPass;

        % Compute the current costs for all destination RF to maintain their current inputs 
        % theCostComponents is an [nDestinationRFs x cost components]  matrix
        theCostComponentsMatrix = obj.totalInputMaintenanceCost();

        % Accumulate cost sequences
        s = theCostComponentsMatrix(1,:);
        idx = find(s > -99);
        netCostSequences = cat(1, netCostSequences, mean(theCostComponentsMatrix(idx,:),1));

        % Visualize convergence
        MosaicConnector.visualizeConvergenceSequence(currentPass, ...
            netCostSequences, obj.costComponentNames(), ...
            netSwaps, obj.wiringParams.maxPassesNum, ...
            'input swaps', 5060);

        % Determine whether convergence was achieved
        convergenceAchieved = MosaicConnector.convergenceAchieved(netCostSequences(:,1));

    end % while

    obj.connectivityMatrix = lastConnectivityMatrix;
    obj.destinationRFspacingsFromCentroids = lastDestinationRFspacingsFromCentroids;

    hFig = figure(5060);
    NicePlot.exportFigToPDF('swapOptimization.pdf', hFig, 300);

    if (generateProgressVideo)
        videoOBJ.close();
    end

    % Save the metaDataStuct for this stage
    if (obj.saveIntermediateConnectivityStagesMetaData)
        obj.updateIntermediateMetaDataStructs();
    end
  
    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.intermediateFigureHandles{numel(obj.intermediateFigureHandles)+1} = ...
            obj.visualizeCurrentConnectivity(1005);
    end

end
