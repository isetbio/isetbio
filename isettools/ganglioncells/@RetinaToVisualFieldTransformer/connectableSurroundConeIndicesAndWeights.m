function [surroundConeIndices, surroundConeWeights] = connectableSurroundConeIndicesAndWeights(...
        surroundConeIndices, surroundConeWeights, modelConstants)

    % Select only the connectable surround cones
    allSurroundConeIndices = surroundConeIndices; 
    allSurroundConeTypes = modelConstants.theConeMosaic.coneTypes(allSurroundConeIndices);
    connectableSurroundConeIndices = [];
    for iConeType = 1:numel(modelConstants.surroundConnectableConeTypes)
        connectableSurroundConeIndices = cat(1, connectableSurroundConeIndices, find(allSurroundConeTypes == modelConstants.surroundConnectableConeTypes(iConeType)));
    end % iConeType
    surroundConeIndices = surroundConeIndices(connectableSurroundConeIndices);
    surroundConeWeights = surroundConeWeights(connectableSurroundConeIndices);
end
