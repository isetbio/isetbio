function updateIntermediateMetaDataStructs(obj,phaseDescriptor, netCostSequences, netTransfers)
    metaDataStruct = {};

    metaDataStruct.phaseDescriptor = phaseDescriptor;
    metaDataStruct.netCostSequences = netCostSequences;
    metaDataStruct.netTransfers = netTransfers;
    metaDataStruct.costComponentNames = obj.costComponentNames();

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
