function [surroundConeIndices, surroundConeWeights, ...
          nonConnectableSurroundConeIndices, ...
          nonConnectableSurroundConeWeights] = connectableSurroundConeIndicesAndWeights(...
                surroundConeIndices, surroundConeWeights, modelConstants)

    % Select only the connectable surround cones
    allSurroundConeIndices = surroundConeIndices; 
    allSurroundConeTypes = modelConstants.inputConeMosaicConeTypes(allSurroundConeIndices);

    connectableSurroundConeIndices = [];
    for iConeType = 1:numel(modelConstants.surroundConnectableConeTypes)
        connectableSurroundConeIndices = cat(1, connectableSurroundConeIndices, find(allSurroundConeTypes == modelConstants.surroundConnectableConeTypes(iConeType)));
    end % iConeType

    % Non-connectable surround cone indices and weights
    nonConnectableSurroundConeIndices = setdiff(1:numel(allSurroundConeTypes), connectableSurroundConeIndices);
    nonConnectableSurroundConeWeights = surroundConeWeights(nonConnectableSurroundConeIndices);

    % Connectable surround cone indices and weights
    surroundConeIndices = surroundConeIndices(connectableSurroundConeIndices);
    surroundConeWeights = surroundConeWeights(connectableSurroundConeIndices);
end