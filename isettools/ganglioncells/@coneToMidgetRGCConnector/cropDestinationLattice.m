function cropDestinationLattice(obj)
    % Crop the destination lattice (RGC positions) so that we have enough space 
    % for the surround cones

    % Compute the max surround radius
    maxSurroundRadiusMicrons = maxSurroundRadiusMicronsAtEccentricity(obj);

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


function surroundRadiusMicrons = maxSurroundRadiusMicronsAtEccentricity(obj)
    
    % Find max ecc degs
    [maxRadialEccMicrons,idx] = max(sqrt(sum(obj.sourceLattice.RFpositionsMicrons.^2,2)));
    eccDegs = obj.sourceLattice.MMsToDegsConversionFunction(1e-3*obj.sourceLattice.RFpositionsMicrons(idx,:));


    % Instantiate RetinaToVisualFieldTrasformer with the Artal database
    xFormer = RetinaToVisualFieldTransformer('ZernikeDataBase', 'Artal2012');

    % Choose a subject with average optics
    subjectRankOrder = 30;
    subjectRankingEye = 'right eye';
    subjID = xFormer.subjectWithRankInEye(subjectRankOrder, subjectRankingEye);

    analyzedEye = 'right eye';
    pupilDiameterMM = 3.0;
    dStruct = xFormer.estimateConeCharacteristicRadiusInVisualSpace(...
                analyzedEye, eccDegs, subjID, pupilDiameterMM, '');

    anatomicalConeCharacteristicRadiusDegs = dStruct.anatomicalConeCharacteristicRadiusDegs;
    visualConeCharacteristicRadiusDegs = dStruct.visualConeCharacteristicRadiusDegs;

    % Estimate surround radius, assuming 1 cone in the RF center and a
    % Rs/Rc ratio = 6;
    RsRcRatio = 6;
    characteristicRadiiLimit = 2;
    surroundRadiusDegs = characteristicRadiiLimit * (visualConeCharacteristicRadiusDegs*RsRcRatio);
    surroundRadiusMicrons = 1e3 * obj.sourceLattice.DegsToMMsConversionFunction(surroundRadiusDegs);
end
