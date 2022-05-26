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
        videoOBJ = VideoWriter('Step5', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    maxPassesNum = obj.wiringParams.maxSwapPassesNum;
    currentPass = 0;
    swapsInCurrentPass = 1;
    
    while (currentPass < maxPassesNum) && (swapsInCurrentPass>0)
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

        hProgressBar = waitbar(0,'Swapping cone phase ...');

        for iRGC = 1:numel(sortedRGCindices)
            theSourceRGCindex = sortedRGCindices(iRGC);

            % Update progress bar
            waitbar(iRGC/numel(sortedRGCindices), hProgressBar, ...
                sprintf('RGC #%d/%d (PASS %d/%d)', iRGC,numel(sortedRGCindices), currentPass, maxPassesNum));

            % Find up to N neigboring RGCs that are no farther than k x RGC separation
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

            for iNearbyRGC = 1:numel(theNeighboringRGCindices)
                theNearbyRGCindex = theNeighboringRGCindices(iNearbyRGC);
                allNeighboringRGCsInputConeIndices{iNearbyRGC} = find(obj.coneConnectivityMatrix(:,theNearbyRGCindex)>0);
                allNeighboringRGCsInputConeWeights{iNearbyRGC} = full(obj.coneConnectivityMatrix(allNeighboringRGCsInputConeIndices{iNearbyRGC},theNearbyRGCindex));
            end % iNearbyRGC

            % See if there is a beneficial swap
            beneficialSwapWasFound = obj.optimizeSwappingOfConeInputs(...
                    theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
                    theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

            if (beneficialSwapWasFound)
                swapsInCurrentPass = swapsInCurrentPass + 1;
                if (generateProgressVideo)
                    % Visualize current connectivity
                    hFig = obj.visualizeCurrentConnectivityState(1005);
                    videoOBJ.writeVideo(getframe(hFig));
                end
            end

        end % iRGC

        if (swapsInCurrentPass>0)
            fprintf('*** There were %d cone swaps during pass %d. Going for another pass.\n', swapsInCurrentPass, currentPass);
        else
            fprintf('No need to do more than #%d passes.\n', currentPass);
        end

        close(hProgressBar);

    end % while loop

    if (generateProgressVideo)
        videoOBJ.close();
    end
end