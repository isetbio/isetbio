function connectSourceRFsToDestinationRFsBasedOnLocalDensities(obj)
% For each destination RF try to connect to N source RFs that are not further than 1
% destination RF separation away and that are not already connected to
% another destination RF. N is the sourceToDestinationDensityRatio
% The obj.connectivityMatrix sparse matrix is initialized here with the connections 
% established at this step.

    % Initialize centroids. No inputs so set them all to inf
    sourceRFsNum = size(obj.sourceLattice.RFpositionsMicrons,1);
    destinationRFsNum = size(obj.destinationLattice.RFpositionsMicrons,1);
   
    obj.destinationRFcentroidsFromInputs = inf(destinationRFsNum,2);

    % Compute the source-to-destination density ratio map at the current
    % RFpos of the destination lattice
    densityRatioMapAtAllDestinationRFpos = obj.sourceToDestinationDensityRatioMap();

    fprintf('Will connect %d input RFs to %d destination RFs with input:destination density ratios in the range [%1.2f - %1.2f]\n', ...
        sourceRFsNum, destinationRFsNum, ...
        min(densityRatioMapAtAllDestinationRFpos(:)), max(densityRatioMapAtAllDestinationRFpos(:)));

    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestDestinationRFIndices = [];
    connectedSourceRFIndices = [];

    [connectedSourceRFIndices, nearestDestinationRFIndices] = doIt(obj, ...
        densityRatioMapAtAllDestinationRFpos, connectedSourceRFIndices, nearestDestinationRFIndices);
    
    % Initialize the [sourceRFsNum x destinationRFsNum] sparse connectivity matrix
    weights = ones([1 numel(connectedSourceRFIndices)]);
    obj.connectivityMatrix = sparse(...
        connectedSourceRFIndices, nearestDestinationRFIndices, weights, ...
        sourceRFsNum, destinationRFsNum);

    % Update the input-based destination RF centroids
    obj.updateDestinationCentroidsFromInputs(unique(nearestDestinationRFIndices));

    % Save the metaDataStuct for this stage
    if (obj.saveIntermediateConnectivityStagesMetaData)
        phaseDescriptor = 'local density based connections';
        obj.updateIntermediateMetaDataStructs(phaseDescriptor, [], []);
    end

    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.intermediateFigureHandles{numel(obj.intermediateFigureHandles)+1} = ...
            obj.visualizeCurrentConnectivity(1001);
    end

end



function [connectedSourceRFIndices, nearestDestinationRFIndices] = doIt(obj, ...
        densityRatioMapAtAllDestinationRFpos, connectedSourceRFIndices, nearestDestinationRFIndices)

    destinationRFsNum = size(obj.destinationLattice.RFpositionsMicrons,1);
    
    for iDestinationRF = 1:destinationRFsNum 

        % Find the local source:destination densityRatioMap  closest source RF to each destination RF
        [distances, idx] = MosaicConnector.pdist2(...
            obj.sourceLattice.RFpositionsMicrons(obj.connectableSourceRFindices,:), ...
            obj.destinationLattice.RFpositionsMicrons(iDestinationRF,:), ...
            'smallest', max([1 floor(densityRatioMapAtAllDestinationRFpos(iDestinationRF))]));

        closestSourceRFIndices = obj.connectableSourceRFindices(idx);
        if (isempty(closestSourceRFIndices))
            continue;
        end

        % Find which source RFs are not more than the local search region from the destination RF.
        % Search region is: MosaicConnector.maxFractionOfRFspacingDuringInitialNeighborSearch x local RF spacing
        withinInitialSearchRegion = distances <= ...
            MosaicConnector.maxFractionOfRFspacingDuringInitialNeighborSearch * ...
            obj.destinationLattice.RFspacingsMicrons(iDestinationRF);

        % Find which of these closest source RFs are not already connected to
        % another destination RF and are within the local search region
        idx = find((~ismember(closestSourceRFIndices, connectedSourceRFIndices)) & (withinInitialSearchRegion));

        if (isempty(idx))
            continue;
        end
        closestSourceRFIndices = closestSourceRFIndices(idx);

        % Accumulate indices for sparse array construction 
        connectedSourceRFIndices = cat(1, connectedSourceRFIndices, closestSourceRFIndices);
        nearestDestinationRFIndices = cat(1, nearestDestinationRFIndices, iDestinationRF*ones(numel(closestSourceRFIndices),1));
    end % iDestinationRF
end