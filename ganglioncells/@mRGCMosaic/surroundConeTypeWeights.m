function [theSurroundConeTypeWeights, theExclusiveSurroundConeTypeWeights, theSurroundConeTypes, theSurroundConeIndices] = ...
    surroundConeTypeWeights(obj, theRGCindex, theCenterConeIndices)

    % Retrieve the cell's surround cone indices and weights
    connectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, theRGCindex)));
    theSurroundConeIndices = find(connectivityVector > 0.0001);

    % Retrieve the cell's surround cone types
    theSurroundConeTypes = obj.inputConeMosaic.coneTypes(theSurroundConeIndices);

    % Compute the net surround weights for each cone type
    theSurroundConeTypeWeights = netSurroundWeightsForEachConeType(...
        connectivityVector, theSurroundConeIndices, theSurroundConeTypes);
    
    % If a non-empty centerConeIndices is passed, also compute the
    % exclusive surround cone type weights
    if (~isempty(theCenterConeIndices))
        [~, idx] = setdiff(theSurroundConeIndices, theCenterConeIndices);
        theExclusiceSurroundConeIndices = theSurroundConeIndices(idx);
        theExclusiceSurroundConeTypes = theSurroundConeTypes(idx);

        % Compute the net surround weights for each cone type
        theExclusiveSurroundConeTypeWeights = netSurroundWeightsForEachConeType(...
            connectivityVector, theExclusiceSurroundConeIndices, theExclusiceSurroundConeTypes);
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
