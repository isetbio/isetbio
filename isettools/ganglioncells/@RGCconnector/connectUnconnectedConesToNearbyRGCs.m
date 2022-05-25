function connectUnconnectedConesToNearbyRGCs(obj, varargin)
% For each unconnected cone try to connect it to the closest RGC 
% that is not further than 1 RGC separation away, and update the
% obj.coneConnectivityMatrix sparse matrix accordingly
    
     % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;

    % Find the indices of cones that are inside the boundary defined by cones 
    % that are already connected to some RGC
    connectedConeIndices = find(squeeze(sum(obj.coneConnectivityMatrix,2)) > 0);
    conesIndicesInsideBoundaryOfAlreadyConnectedCones = ...
        RGCconnector.pointsInsideBoundaryDefinedBySelectedPoints(...
            obj.inputConeMosaic.coneRFpositionsMicrons, connectedConeIndices);


    % Find indices of LM cones that are not already connected to some RGC
    % and are inside the boundary defined by cones that are connected to
    % some RGC
    unconnectedLMconeIndices = find(...
        (conesIndicesInsideBoundaryOfAlreadyConnectedCones) & ...
        (squeeze(sum(obj.coneConnectivityMatrix,2)) == 0) & ...
        (obj.inputConeMosaic.coneTypes ~= cMosaic.SCONE_ID));
    unconnectedLMconesNum = numel(unconnectedLMconeIndices);


    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step2', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    for iCone = 1:numel(unconnectedLMconeIndices)
        theConeIndex = unconnectedLMconeIndices(iCone);

        % Find the closest RGC to each this cone
        [distance, theTargetRGCindex] = RGCconnector.pdist2(...
            obj.RGCRFcentroidsFromInputs, ...
            obj.inputConeMosaic.coneRFpositionsMicrons(theConeIndex,:), ...
            '', ...
            'smallest', 1);

        normDistance = distance/obj.RGCRFspacingsMicrons(theTargetRGCindex);

        if (normDistance > 1.5)
            error('WIll not connect cone %d to any RGC. Closest RGC has a distance of %2.2f\n', iCone, normDistance);
     
        else
            %fprintf('Connecting %d to %d RGC which is %2.2f distance away\n', iCone,theTargetRGCindex, normDistance);
            
            % Update connectivity matrix
            obj.coneConnectivityMatrix(theConeIndex, theTargetRGCindex) = 1;

            % Update the centroid
            obj.updateCentroidsFromInputs(theTargetRGCindex);

            if (generateProgressVideo)
                % Visualize current connectivity
                hFig = obj.visualizeCurrentConnectivityState(1001);
                videoOBJ.writeVideo(getframe(hFig));
            end
        end
    end % iCone

    % Update spacings based on new centroids
    obj.updateRGCRFspacingsBasedOnCurrentCentroids();

    if (generateProgressVideo)
        videoOBJ.close();
    end

end