function transferConesBetweenNearbyRGCsWithUnbalancedInputNumerosities(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('optimizationCenter', 'patchCenter', @(x)(ismember(x, {'patchCenter', 'visualFieldCenter'})));

    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;
    optimizationCenter = p.Results.optimizationCenter;

    minConeInputsNum = 3;
    maxConeInputsNum = 5;

    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step3', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    maxPassesNum = 10;
    currentPass = 0;
    tranfersInCurrentPass = 1;
    
    while (currentPass < maxPassesNum) && (tranfersInCurrentPass>0)

        % Update pass number
        currentPass = currentPass + 1;
        tranfersInCurrentPass = 0;

        % We only consider RGCs that have between [minConeInputsNum maxConeInputsNum] input cones
        ss = squeeze(sum(obj.coneConnectivityMatrix,1));
        targetedRGCindices = find((ss >= minConeInputsNum) & (ss <= maxConeInputsNum));
        if (isempty(targetedRGCindices))
            continue;
        end

        % Retrieve their centroids
        targetedRGCCentroids = obj.RGCRFcentroidsFromInputs(targetedRGCindices,:);
    
        % Sort RGCs according to their ecc
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
            % Find the indices of cone inputs to this RGC
            theSourceRGCindex = sortedRGCindices(iRGC);
            theSourceRGCinputConeIndices = find(squeeze(obj.coneConnectivityMatrix(:,theSourceRGCindex))>0);
            theSourceRGCinputConesNum = numel(theSourceRGCinputConeIndices);
            theSourceRGCinputConeWeights = full(obj.coneConnectivityMatrix(theSourceRGCinputConeIndices,theSourceRGCindex));
    
            % Find up to N neigboring RGCs that are no farther than k x RGC separation
            nearbyRGCindices = obj.neihboringRGCindices(theSourceRGCindex);
    
            if (~isempty(nearbyRGCindices))
                % Find the difference in # of input cones between source RGC and neigboring RGCs
                inputConeIndicesAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
                inputConeWeightsAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
                inputConeTypesAllNeighboringRGCs = cell(1,numel(nearbyRGCindices));
    
                differenceInConeInputs = zeros(1,numel(nearbyRGCindices));
                for iNearbyRGC = 1:numel(nearbyRGCindices)
                    theNearbyRGCindex = nearbyRGCindices(iNearbyRGC);
                    inputConeIndicesAllNeighboringRGCs{iNearbyRGC} = find(obj.coneConnectivityMatrix(:,theNearbyRGCindex)>0);
                    inputConeWeightsAllNeighboringRGCs{iNearbyRGC} = full(obj.coneConnectivityMatrix(inputConeIndicesAllNeighboringRGCs{iNearbyRGC},theNearbyRGCindex));
                    inputConeTypesAllNeighboringRGCs{iNearbyRGC} = obj.inputConeMosaic.coneTypes(inputConeIndicesAllNeighboringRGCs{iNearbyRGC});
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
    
                end
            end % (~isempty(nearbyRGCindices))
        end % iRGC

        if (tranfersInCurrentPass>0)
            fprintf('*** There were %d cone transfers during pass %d. Going for another pass.\n', tranfersInCurrentPass, currentPass);
        else
            fprintf('No need to do more than #%d passes.\n', currentPass);
        end

    end % pass

    if (generateProgressVideo)
        videoOBJ.close();
    end

end