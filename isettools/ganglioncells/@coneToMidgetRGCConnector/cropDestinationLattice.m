function cropDestinationLattice(obj)
    % Crop the destination lattice (RGC positions) so that we have enough space 
    % for the surround cones

    fprintf(2, '********************************\n')
    fprintf(2,'Need to compute max SurroundDiameterMicrons for mRGCs at this ecc\n')
    fprintf(2, 'Arbitrarily setting this to 20 microns\n');
    fprintf(2, '********************************\n')

    maxSurroundDiameterMicrons = 20;
    maxSurroundRadiusMicrons = 0.5*maxSurroundDiameterMicrons;

    minConePosXY = min(obj.sourceLattice.RFpositionsMicrons,[],1);
    maxConePosXY = max(obj.sourceLattice.RFpositionsMicrons,[],1);

    idx = find(...
        (obj.destinationLattice.RFpositionsMicrons(:,1) >= minConePosXY(1) + maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,1) <= maxConePosXY(1) - maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,2) >= minConePosXY(2) + maxSurroundRadiusMicrons) & ...
        (obj.destinationLattice.RFpositionsMicrons(:,2) <= maxConePosXY(2) - maxSurroundRadiusMicrons));


    % Crop RGC positions
    obj.destinationLattice.RFpositionsMicrons = obj.destinationLattice.RFpositionsMicrons(idx,:);
       
    % Update RGC spacings & nearby RGCindices
    [obj.destinationLattice.RFspacingsMicrons, ...
     obj.destinationLattice.nearbyRFindices] = RGCmodels.Watson.convert.positionsToSpacings(obj.destinationLattice.RFpositionsMicrons);

end
