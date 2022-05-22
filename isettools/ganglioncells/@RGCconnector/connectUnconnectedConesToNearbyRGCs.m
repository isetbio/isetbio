function connectUnconnectedConesToNearbyRGCs(obj)
% For each unconnected cone try to connect it to the closest RGC 
% that is not further than 1 RGC separation away, and update the
% obj.coneConnectivityMatrix sparse matrix accordingly
    
    % Find all 0-, 1-, and 2-input RGCs
    zeroInputRGCsNum = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 0));
    oneInputRGCsNum = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 1));
    twoInputRGCsNum = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 2));
    
    unconnectedLMconeIndices = find(...
        (squeeze(sum(obj.coneConnectivityMatrix,2)) == 0) & ...
        (obj.inputConeMosaic.coneTypes ~= cMosaic.SCONE_ID));
    unconnectedLMconesNum = numel(unconnectedLMconeIndices);

    fprintf('There are %d L/M cones that are not connected to any RGC.\n', unconnectedLMconesNum);
    %fprintf('There are %d RGCS that have 0 cone inputs.\n', zeroInputRGCsNum);
    %fprintf('There are %d RGCS that have 1 cone input.\n', oneInputRGCsNum);
    %fprintf('There are %d RGCS that have 2 cone inputs.\n', twoInputRGCsNum);
    
      
    % Video setup
    videoOBJ = VideoWriter('Step2', 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    for iCone = 1:numel(unconnectedLMconeIndices)
        theConeIndex = unconnectedLMconeIndices(iCone);

        % Find the closest RGC to each this cone
        [distance, theTargetRGCindex] = RGCconnector.pdist2(...
            obj.RGCRFcentroidsFromInputs, ...
            obj.inputConeMosaic.coneRFpositionsMicrons(theConeIndex,:), ...
            '', ...
            'smallest', 1);

        normDistance = distance/obj.RGCRFspacingsMicrons(theTargetRGCindex);
        if (normDistance > 1)
            fprintf(2,'WIll not connect cone %d to any RGC. Closest RGC has a distance of %2.2f\n', iCone, normDistance);
           
        else
            fprintf('Connecting %d to %d RGC which is %2.2f distance away\n', iCone,theTargetRGCindex, normDistance);
            
            % Update connectivity matrix
            obj.coneConnectivityMatrix(theConeIndex, theTargetRGCindex) = 1;

            % Update the centroid
            obj.updateCentroidsFromInputs(theTargetRGCindex);

            % Visualize current connectivity
            hFig = obj.visualizeCurrentConnectivityState(1001);
            videoOBJ.writeVideo(getframe(hFig));
        end

    end % iCone
    videoOBJ.close()


end