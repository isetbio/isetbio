function debug_overwritePositionConnectivityData(obj, theStageMetaDataStruct)

        obj.rgcRFcenterConeConnectivityMatrix = theStageMetaDataStruct.connectivityMatrix;
    
        % Overwrite positions/spacings in microns
        obj.rgcRFpositionsMicrons = theStageMetaDataStruct.destinationLattice.RFpositionsMicrons;
        obj.rgcRFspacingsMicrons = reshape(theStageMetaDataStruct.destinationLattice.RFspacingsMicrons, [1 numel(theStageMetaDataStruct.destinationLattice.RFspacingsMicrons)]);
    
        % And in degrees
        obj.rgcRFpositionsDegs = obj.inputConeMosaic.customMMsToDegsConversionFunction(obj.rgcRFpositionsMicrons/1000);
        
        % Set the rf spacings in degrees
        RFeccentricityMicrons = sqrt(sum((obj.rgcRFpositionsMicrons).^2,2));
        obj.rgcRFspacingsDegs = reshape(obj.inputConeMosaic.sizeMicronsToSizeDegreesForCmosaic(obj.rgcRFspacingsMicrons, reshape(RFeccentricityMicrons, size(obj.rgcRFspacingsMicrons))), ...
                                    [1 numel(theStageMetaDataStruct.destinationLattice.RFspacingsMicrons)]);
    
        % Reset the visualization cache
        obj.visualizationCache.rfCenterPatchData = [];
    end