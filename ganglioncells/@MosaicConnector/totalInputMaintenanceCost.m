function theCostComponentsMatrix = totalInputMaintenanceCost(obj)

    destinationRFsNum = size(obj.connectivityMatrix,2);

    % Initialize the theCostComponentsMatrix
    iDestinationRF = 1;
    inputIndices = find(squeeze(obj.connectivityMatrix(:, iDestinationRF))>0);
    inputWeights = full(obj.connectivityMatrix(inputIndices,iDestinationRF));
    destinationRFspacing = obj.destinationRFspacingsFromCentroids(iDestinationRF);

    tmp = obj.inputMaintenanceCost(inputIndices, inputWeights, destinationRFspacing);
    theCostComponentsMatrix = zeros(destinationRFsNum, numel(tmp));

    cm = obj.connectivityMatrix;
    
    % Compute current costs
    for iDestinationRF = 1:destinationRFsNum

        inputIndices = find(squeeze(cm(:, iDestinationRF))>0);
        inputWeights = full(cm(inputIndices,iDestinationRF));
        destinationRFspacing = obj.destinationRFspacingsFromCentroids(iDestinationRF);

        theCostComponentsMatrix(iDestinationRF, :) = ...
            obj.inputMaintenanceCost(inputIndices, inputWeights, destinationRFspacing);
    end
end

