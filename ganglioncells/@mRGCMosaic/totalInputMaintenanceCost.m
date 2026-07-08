function theCostComponentsMatrix = totalInputMaintenanceCost(obj, ...
    chromaticSpatialVarianceTradeoff, spatialVarianceMetric)

    % Initialize the theCostComponentsMatrix
    theCostComponentsMatrix = zeros(obj.rgcsNum, 3);

    % Retrieve the cone pooling matrix
    if (~isempty(obj.rgcRFcenterConeConnectivityMatrix))
        cm = obj.rgcRFcenterConeConnectivityMatrix;
    else
        cm = obj.rgcRFcenterConePoolingMatrix;
    end

    % Compute costs for each RGC
    for iRGC = 1:obj.rgcsNum
        inputIndices = find(squeeze(cm(:, iRGC))>0);
        inputWeights = full(cm(inputIndices,iRGC));
        rgcRFspacing = obj.rgcRFspacingsMicrons(iRGC);

        theCostComponentsMatrix(iRGC, :) = inputMaintenanceCost(...
            obj.inputConeMosaic, ...
            inputIndices, inputWeights, rgcRFspacing, ...
            chromaticSpatialVarianceTradeoff, spatialVarianceMetric);
    end

end



function theCostComponents = inputMaintenanceCost(...
    inputConeMosaic, ...
    inputConeIndices, inputConeWeights, rgcRFspacing, ...
    chromaticSpatialVarianceTradeoff, spatialVarianceMetric)

    if (isempty(inputConeIndices))
        %fprintf(2,'\t\t* * * * - - - > Cost to maintain 0 cones: -99\n');
        theCostComponents = [-99 -99 -99];
        return;
    end

    if (numel(inputConeIndices) == 1)
        % zero cost to maintain 1 cone
        theCostComponents = [0 0 0];
        return;
    else
        % Retrieve the input cone positions
        inputConePositions = inputConeMosaic.coneRFpositionsMicrons(inputConeIndices,:);

        % Retrieve the input cone types
        inputConeTypes = inputConeMosaic.coneTypes(inputConeIndices);

        % Spatial variance cost
        spatialVarianceCost = coneToMidgetRGCConnector.spatialVarianceCost(...
            spatialVarianceMetric, inputConeWeights, inputConePositions, rgcRFspacing);

        % Chromatic variance cost
        chromaticVarianceCost = coneToMidgetRGCConnector.chromaticVarianceCost(inputConeWeights, inputConeTypes);

        % Total cost
        w = chromaticSpatialVarianceTradeoff;
        cost = w * spatialVarianceCost + (1-w)*chromaticVarianceCost;
    end

    % Form vector with cost components
    theCostComponents = [cost spatialVarianceCost chromaticVarianceCost];

end
