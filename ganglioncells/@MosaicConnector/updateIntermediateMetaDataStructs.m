function updateIntermediateMetaDataStructs(obj)
    metaDataStruct = {};
    metaDataStruct.connectivityMatrix = obj.connectivityMatrix;
    metaDataStruct.sourceLattice.RFpositionsMicrons = obj.sourceLattice.RFpositionsMicrons;
    metaDataStruct.sourceLattice.RFspacingsMicrons = obj.sourceLattice.RFspacingsMicrons;

    if (isempty(obj.destinationRFspacingsFromCentroids))
        metaDataStruct.destinationLattice.RFpositionsMicrons = obj.destinationLattice.RFpositionsMicrons;
        metaDataStruct.destinationLattice.RFspacingsMicrons = obj.destinationLattice.RFspacingsMicrons;
    else
        metaDataStruct.destinationLattice.RFpositionsMicrons = obj.destinationRFcentroidsFromInputs;
        metaDataStruct.destinationLattice.RFspacingsMicrons = obj.destinationRFspacingsFromCentroids;
    end
    obj.intermediateMetaDataStructs{numel(obj.intermediateMetaDataStructs)+1} = metaDataStruct;
end
