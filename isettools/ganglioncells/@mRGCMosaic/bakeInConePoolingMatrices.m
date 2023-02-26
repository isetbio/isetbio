function bakeInConePoolingMatrices(obj, centerConePoolingMatrix, surroundConePoolingMatrix)

    % The centerConeConnectivityMatrix is no longer operational
    % Empty it to save space
    obj.rgcRFcenterConeConnectivityMatrix = [];

    obj.rgcRFcenterConePoolingMatrix = centerConePoolingMatrix;
    obj.rgcRFsurroundConePoolingMatrix = surroundConePoolingMatrix;
end
