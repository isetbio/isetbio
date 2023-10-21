function [theSurroundConeTypeWeights, theExclusiveSurroundConeTypeWeights] = ...
    surroundConeTypeWeights(obj, theRGCindex, theCenterConeIndices)

    % Retrieve the cell's surround cone indices and weights
    connectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, theRGCindex)));
    surroundConeIndices = find(connectivityVector > 0.0001);

    % Retrieve the cell's surround cone types
    surroundConeTypes = obj.inputConeMosaic.coneTypes(surroundConeIndices);

    % Compute the net surround weights for each cone type
    theSurroundConeTypeWeights = netSurroundWeightsForEachConeType(...
        connectivityVector, surroundConeIndices, surroundConeTypes);
    
    % If a non-empty centerConeIndices is passed, also compute the
    % exclusive surround cone type weights
    if (~isempty(theCenterConeIndices))
        [~, idx] = setdiff(surroundConeIndices, theCenterConeIndices);
        surroundConeIndices = surroundConeIndices(idx);
        surroundConeTypes = surroundConeTypes(idx);

        % Compute the net surround weights for each cone type
        theExclusiveSurroundConeTypeWeights = netSurroundWeightsForEachConeType(...
            connectivityVector, surroundConeIndices, surroundConeTypes);
    else
        theExclusiveSurroundConeTypeWeights = [];
    end

end

function theSurroundConeTypeWeights = netSurroundWeightsForEachConeType(connectivityVector, surroundConeIndices, surroundConeTypes)
    % Retrieve all surround weights
    surroundConeWeights = connectivityVector(surroundConeIndices);

    idx = find(surroundConeTypes == cMosaic.LCONE_ID);
    netSurroundLconeWeight = sum(surroundConeWeights(idx));

    idx = find(surroundConeTypes == cMosaic.MCONE_ID);
    netSurroundMconeWeight = sum(surroundConeWeights(idx));

    idx = find(surroundConeTypes == cMosaic.SCONE_ID);
    netSurroundSconeWeight = sum(surroundConeWeights(idx));

    theSurroundConeTypeWeights(cMosaic.LCONE_ID) = netSurroundLconeWeight;
    theSurroundConeTypeWeights(cMosaic.MCONE_ID) = netSurroundMconeWeight;
    theSurroundConeTypeWeights(cMosaic.SCONE_ID) = netSurroundSconeWeight;
end
