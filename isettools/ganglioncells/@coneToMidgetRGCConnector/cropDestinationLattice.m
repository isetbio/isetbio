function cropDestinationLattice(obj)
    % Crop the destination lattice (RGC positions) so that we have enough space 
    % for the surround cones

    % Compute the max surround radius
    maxSurroundRadiusMicrons = obj.sourceLattice.metaData.midgetRGCSurroundRadiusMicronsAtMaxEccentricityGivenOptics;
    minConePosXY = min(obj.sourceLattice.RFpositionsMicrons,[],1);
    maxConePosXY = max(obj.sourceLattice.RFpositionsMicrons,[],1);

    idx = find(...
        (obj.destinationLattice.RFpositionsMicrons(:,1) >= minConePosXY(1) + maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,1) <= maxConePosXY(1) - maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,2) >= minConePosXY(2) + maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,2) <= maxConePosXY(2) - maxSurroundRadiusMicrons));


    if (isempty(idx))
        error('Zero destination RFs left after cropping the RGC lattice to account for the surround radius. Increase size.')
    end

    % Crop RGC positions
    obj.destinationLattice.RFpositionsMicrons = obj.destinationLattice.RFpositionsMicrons(idx,:);
    

    % Update RGC spacings & nearby RGCindices
    obj.destinationLattice.RFspacingsMicrons = obj.destinationLattice.RFspacingsMicrons(idx);
    obj.destinationLattice.nearbyRFindices = obj.destinationLattice.nearbyRFindices(idx,:);
end