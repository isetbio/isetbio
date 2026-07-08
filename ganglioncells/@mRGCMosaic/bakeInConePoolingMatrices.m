function bakeInConePoolingMatrices(obj, centerConePoolingMatrix, surroundConePoolingMatrix)

    % Sanity checks
    assert(all(size(centerConePoolingMatrix) == size(obj.rgcRFcenterConeConnectivityMatrix)), ...
        sprintf('Mismatch in size of rgcRFcenterConePoolingMatrix'));

    assert(all(size(surroundConePoolingMatrix) == size(obj.rgcRFcenterConeConnectivityMatrix)), ...
        sprintf('Mismatch in size of rgcRFsurroundConePoolingMatrix'));

    % Remove intermediate connectivity data 
    obj.intermediateMetaDataStructs = {};
    obj.mosaicConnectorIntermediateFigureHandles = {};

    % The centerConeConnectivityMatrix is no longer operational
    % Empty it to save space
    obj.rgcRFcenterConeConnectivityMatrix = [];

    obj.rgcRFcenterConePoolingMatrix = centerConePoolingMatrix;
    obj.rgcRFsurroundConePoolingMatrix = surroundConePoolingMatrix;
end
