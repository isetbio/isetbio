function connectUnconnectedSourceRFsToClosestDestinationRF(obj)

    % Find the indices of source lattice RFs that are inside the boundary defined by
    % source RFs that have been connected to some destination RF
    connectedSourceRFIndices = find(squeeze(sum(obj.connectivityMatrix,2)) > 0);

    sourceRFIndicesInsideBoundaryOfConnectedSourceRFs = pointsInsideBoundaryDefinedBySelectedPoints(...
            obj.sourceLattice.RFpositionsMicrons, connectedSourceRFIndices);

    
    % Find indices of connectable source RFs that are not already connected
    % and are inside the above boundary, i.e. the intersection of
    % sourceRFIndicesInsideBoundaryOfConnectedSourceRFs
    % obj.connectableSourceRFindices

    idx = find(squeeze(sum(obj.connectivityMatrix(obj.connectableSourceRFindices,:),2)) == 0);
    currentlyUnconnectedButConnectableSourceRFindices = obj.connectableSourceRFindices(idx);

    sourceRFindicesToBeConnected = intersect(...
        currentlyUnconnectedButConnectableSourceRFindices, ...
        sourceRFIndicesInsideBoundaryOfConnectedSourceRFs);
    
    for iSourceRF = 1:numel(sourceRFindicesToBeConnected)
        theSourceRFindex = sourceRFindicesToBeConnected(iSourceRF);

        % Find the closest RGC to each this cone
        [distance, theTargetDestinationRFindex] = MosaicConnector.pdist2(...
            obj.destinationRFcentroidsFromInputs, ...
            obj.sourceLattice.RFpositionsMicrons(theSourceRFindex,:), ...
            'smallest', 1);

        normDistance = distance/obj.destinationLattice.RFspacingsMicrons(theTargetDestinationRFindex);

        if (normDistance > obj.wiringParams.maxNeighborNormDistance)
            error('WIll not connect sourceRF %d to any destination RF. Closest destination RF has a distance of %2.2f\n', theSourceRFindex, normDistance);
        else
            %fprintf('Connecting sourceRF %d to destination RF %d which is %2.2f distance away\n', ...
            %    theSourceRFindex, theTargetDestinationRFindex, normDistance);
 
            % Update connectivity matrix
            obj.connectivityMatrix(theSourceRFindex, theTargetDestinationRFindex) = 1;

            % Update the centroid
            obj.updateDestinationCentroidsFromInputs(theTargetDestinationRFindex);
        end

    end % iSourceRF

    % Update the centroid-based RF spacings since the centroids have changed
    obj.updateDestinationRFspacingsBasedOnCentroids();

    % Save the metaDataStuct for this stage
    if (obj.saveIntermediateConnectivityStagesMetaData)
        phaseDescriptor = 'connecting unconnected source RFs to nearest destination RF';
        obj.updateIntermediateMetaDataStructs(phaseDescriptor, [], []);
    end

    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.intermediateFigureHandles{numel(obj.intermediateFigureHandles)+1} = ...
            obj.visualizeCurrentConnectivity(1002);
    end
end

function [insideBoundaryPointIndices, onBoundaryPointIndices] = pointsInsideBoundaryDefinedBySelectedPoints(...
    allPointPositions, selectedPointIndices)

    % Find indices of points are inside the boundary defined by a select subset of points

    allXcoords = squeeze(allPointPositions(:,1));
    allYcoords = squeeze(allPointPositions(:,2));
    
    shrinkFactor = 1.0;
    idx = boundary(allXcoords(selectedPointIndices), allYcoords(selectedPointIndices), shrinkFactor);

    boundingPolygonXcoords = allXcoords(selectedPointIndices(idx));
    boundingPolygonYcoords = allYcoords(selectedPointIndices(idx));
    [in,on] = inpolygon(allXcoords, allYcoords, boundingPolygonXcoords, boundingPolygonYcoords);

    insideBoundaryPointIndices = find((in == 1) | (on == 1));
    onBoundaryPointIndices = find(on == 1);
end
