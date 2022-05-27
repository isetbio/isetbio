function [totalCost, spatialCost, chromaticCost] = computeInputMaintenanceCostAcrossEntireMosaic(obj)
    rgcsNum = size(obj.coneConnectivityMatrix,2);
    totalCost = zeros(1, rgcsNum);
    spatialCost = zeros(1, rgcsNum);
    chromaticCost = zeros(1, rgcsNum);
    
    if (isempty(obj.localRGCRFspacingsMicrons))
        error('local RGCRFspacings microns must be set by now')
    else
        RGCRFspacings = obj.localRGCRFspacingsMicrons;
    end

    % Compute current costs
    parfor iRGC = 1:rgcsNum
        % Indices of input cones
        connectedConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
        % Weights of input cones
        inputConeWeights = full(obj.coneConnectivityMatrix(connectedConeIndices, iRGC));
        % Compute the costs
        [totalCost(iRGC), spatialCost(iRGC), chromaticCost(iRGC)] = obj.costToMaintainInputs(...
            connectedConeIndices, inputConeWeights, ...
            RGCRFspacings(iRGC));
    end
end