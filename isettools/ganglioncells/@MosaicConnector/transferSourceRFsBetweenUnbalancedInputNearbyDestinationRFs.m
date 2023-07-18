function transferSourceRFsBetweenUnbalancedInputNearbyDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;

    % We only transfer inputs from destinationRF-A to a nearby destination RF-B
    % if the destinationRF-A has between 3 (MosaicConnector.minConeInputsPerRGCToConsiderTransferToNearbyRGCs)
    % and obj.wiringParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs inputs
    % If it has less than 3 inputs, i.e. 2, then this will be dealt by the
    % swapSourceRFsBetweenNearbyDestinationRFs() method which is called later on
    minInputsNum = MosaicConnector.minSourceInputsToConsiderTransferToNearbyDestinationRF; 
    maxInputsNum = obj.wiringParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs;

    fprintf('Will optimize homogeneity of cone inputs to nearby RGCs for RGCs with center cones num in the range [%d .. %d]\n', ...
        minInputsNum, maxInputsNum);
    ss = squeeze(sum(obj.connectivityMatrix,1));
    targetedDestinationRFindices = find((ss >= minInputsNum) & (ss <= maxInputsNum));
    if (isempty(targetedDestinationRFindices))
        fprintf('There are no destination RFs with 3-6 inputs. This phase will not execute. \n.')
        return;
    end


    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step3', 'MPEG-4');
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

    currentPass = 0;  tranfersInCurrentPass = 1; convergenceAchieved = false;
    while (currentPass < obj.wiringParams.maxPassesNum) && ...
          (tranfersInCurrentPass>0) && ...
          (~convergenceAchieved)

        % Update pass number
        currentPass = currentPass + 1;
        tranfersInCurrentPass = 0;

        % We only consider destination RFs that have between [minInputsNum maxInputsNum] input cones
        ss = squeeze(sum(obj.connectivityMatrix,1));
        targetedDestinationRFindices = find((ss >= minInputsNum) & (ss <= maxInputsNum));
        if (isempty(targetedDestinationRFindices))
            continue;
        end

        % Sort targetedDestinationRF indices (based on their eccentricity &
        % optimization center)
        idx = obj.sortDestinationRFsBasedOnOptimizationCenter(targetedDestinationRFindices);
        sortedDestinationRFindices = targetedDestinationRFindices(idx);

        for iDestinationRF = 1:numel(sortedDestinationRFindices)
            theSourceDestinationRFindex = sortedDestinationRFindices(iDestinationRF);
            
            % Find the indicies of the neigboring destinationRFs
            nearbyDestinationRFindices = obj.indicesOfNeighboringDestinationRFs(theSourceDestinationRFindex);
            if (isempty(nearbyDestinationRFindices))
                continue;
            end

            % Find the indices and weights of cone inputs to this destination RF
            theSourceDestinationRFinputIndices = find(squeeze(obj.connectivityMatrix(:,theSourceDestinationRFindex))>0);
            theSourceDestinationRFinputConesNum = numel(theSourceDestinationRFinputIndices);
            theSourceDestinationRFinputWeights = full(obj.connectivityMatrix(theSourceDestinationRFinputIndices,theSourceDestinationRFindex));
    
            % Find the difference in # of inputs between source
            % destinationRF and nearby destinationRFs
            inputIndicesAllNearbyDestinationRFs = cell(1,numel(nearbyDestinationRFindices));
            inputWeightsAllNearbyDestinationRFs = cell(1,numel(nearbyDestinationRFindices));
            differenceInInputs = zeros(1,numel(nearbyDestinationRFindices));
            for iNearbyDestinationRF = 1:numel(nearbyDestinationRFindices)
                theNearbyDestinationRFindex = nearbyDestinationRFindices(iNearbyDestinationRF);
                inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF} = find(obj.connectivityMatrix(:,theNearbyDestinationRFindex)>0);
                inputWeightsAllNearbyDestinationRFs{iNearbyDestinationRF} = full(obj.connectivityMatrix(inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF},theNearbyDestinationRFindex));
                differenceInInputs(iNearbyDestinationRF) = theSourceDestinationRFinputConesNum - numel(inputIndicesAllNearbyDestinationRFs{iNearbyDestinationRF});
            end % iNearbyDestinationRF

            % Find nearby destinationRFs whose # of inputs is at least 2 less
            maxDiff = max(differenceInInputs);
            if (maxDiff >= 2)
                tranfersInCurrentPass = tranfersInCurrentPass + 1;
                idx = find(differenceInInputs == maxDiff);
                theNearbyDestinationRFInputIndices = inputIndicesAllNearbyDestinationRFs(idx);
                theNearbyDestinationRFInputWeights = inputWeightsAllNearbyDestinationRFs(idx);
                theNearbyDestinationRFindex = nearbyDestinationRFindices(idx);
                obj.optimizeTransferOfInputRFs(...
                    theSourceDestinationRFindex, theSourceDestinationRFinputIndices, theSourceDestinationRFinputWeights, ...
                    theNearbyDestinationRFindex, theNearbyDestinationRFInputIndices, theNearbyDestinationRFInputWeights);

                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivity(1003);
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end %if (maxDiff >= 2)
        end % for iDestinationRF

        % Update the destinationRF spacings based on the updated connectivity
        obj.updateDestinationRFspacingsBasedOnCentroids();

        % Update transfers sequence
        netTransfers(currentPass) = tranfersInCurrentPass;

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
            netTransfers, obj.wiringParams.maxPassesNum, ...
            'input transfers', 5050);

        % Determine whether convergence was achieved
        convergenceAchieved = MosaicConnector.convergenceAchieved(netCostSequences(:,1));

    end % while loop

    hFig = figure(5050);
    NicePlot.exportFigToPDF('transferOptimization.pdf', hFig, 300);

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
            obj.visualizeCurrentConnectivity(1003);
    end
end