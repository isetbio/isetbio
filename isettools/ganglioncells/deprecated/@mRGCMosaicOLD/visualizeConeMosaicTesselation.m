function visualizeConeMosaicTesselation(obj, domain)

    switch (domain)
        case 'microns'
            rgcRFpositions = obj.rgcRFpositionsMicrons;
            rgcRFspacings = obj.rgcRFspacingsMicrons;
            coneRFpositions = obj.inputConeMosaicMetaData.conePositionsMicrons;
            coneRFspacings = obj.inputConeMosaicMetaData.coneSpacingsMicrons;
        case 'degrees'
            rgcRFpositions = obj.rgcRFpositionsDegs;
            rgcRFspacings = obj.rgcRFspacingsDegs;
            coneRFpositions = obj.inputConeMosaicMetaData.conePositionsDegs;
            coneRFspacings = obj.inputConeMosaicMetaData.coneSpacingsDegs;
        otherwise
            error('Unknown visualization domain: ''%s''.', domain);
    end
    
    coneTypes = obj.inputConeMosaicMetaData.coneTypes;
    mRGCmosaic.renderTesselationMap([], [], coneRFpositions, coneRFspacings, coneTypes, ...
        rgcRFpositions, rgcRFspacings, obj.coneConnectivityMatrix, domain);
    
end

