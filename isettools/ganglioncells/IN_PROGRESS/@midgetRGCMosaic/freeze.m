function freeze(obj)
    
    % Set the isFrozen flag
    obj.isFrozen = true;

    % Remove large chunks of data that are of no use for the frozen mosaic
    for iObj = 1:numel(obj.theRTVFobjList)
        obj.theRTVFobjList{iObj}.freeze();
    end

    % Remove the rgcRFcenterConeConnectivityMatrix = [];
    obj.rgcRFcenterConeConnectivityMatrix = [];

    % Remove theMosaicConnectorOBJ
    obj.theMosaicConnectorOBJ = [];
end
