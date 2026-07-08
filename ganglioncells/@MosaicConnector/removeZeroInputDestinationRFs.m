function removeZeroInputDestinationRFs(obj, indicesOfZeroInputDestinationRFs)

    % Crop RGC positions
    allDestinationRFindices = 1:size(obj.destinationLattice.RFpositionsMicrons,1);
    idx = setdiff(allDestinationRFindices, indicesOfZeroInputDestinationRFs);

    % Remove from destinationLattice
    obj.destinationLattice.RFpositionsMicrons = obj.destinationLattice.RFpositionsMicrons(idx,:);
    
    % Update RGC spacings & nearby RGCindices
    obj.destinationLattice.RFspacingsMicrons = obj.destinationLattice.RFspacingsMicrons(idx);
    obj.destinationLattice.nearbyRFindices = obj.destinationLattice.nearbyRFindices(idx,:);

    % Update connectivity matrix
    obj.connectivityMatrix = obj.connectivityMatrix(:,idx);

    % Update destinationRFcentroidsFromInputs
    obj.destinationRFcentroidsFromInputs = obj.destinationRFcentroidsFromInputs(idx,:);

    % Update destinationRFspacingsFromCentroids
    obj.destinationRFspacingsFromCentroids = obj.destinationRFspacingsFromCentroids(idx);
end