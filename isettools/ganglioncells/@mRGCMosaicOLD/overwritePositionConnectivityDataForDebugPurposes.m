function overwritePositionConnectivityDataForDebugPurposes(obj, theStageMetaDataStruct)
    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        error('Overwriting position/connectivity data for an mRGCmosaic at this stage is not allowed.')
    end

    obj.rgcRFcenterConeConnectivityMatrix = theStageMetaDataStruct.connectivityMatrix;

    % Overwrite positions/spacings in microns
    obj.rgcRFpositionsMicrons = theStageMetaDataStruct.destinationLattice.RFpositionsMicrons;
    obj.rgcRFspacingsMicrons = theStageMetaDataStruct.destinationLattice.RFspacingsMicrons;

    % And in degrees
    obj.rgcRFpositionsDegs = obj.inputConeMosaic.customMMsToDegsConversionFunction(obj.rgcRFpositionsMicrons/1000);
    
    % Set the rf spacings in degrees
    RFeccentricityMicrons = sqrt(sum((obj.rgcRFpositionsMicrons).^2,2));
    obj.rgcRFspacingsDegs = obj.inputConeMosaic.sizeMicronsToSizeDegreesForCmosaic(obj.rgcRFspacingsMicrons, reshape(RFeccentricityMicrons, size(obj.rgcRFspacingsMicrons)));


    % Reset the visualization cache
    obj.visualizationCache.rfCenterPatchData = [];

end
