function [coneWeights, rgcPositionsDegsFromConnectivity, synthesizedRFparams] = computeConeWeights(...
    conePositionsDegs,  coneTypes, connectivityMatrix)
% Method to compute weights of connections b/n cones and  center-surround RF subregions
        
    conesNum = size(connectivityMatrix,1);
    rgcsNum = size(connectivityMatrix,2);
    rgcCenterConeInputsNum = zeros(rgcsNum,1);
    rgcPositionsDegsFromConnectivity = zeros(rgcsNum,2);
    
    % Compute number of cone inputs to the RGC RF centers
    for RGCindex = 1:rgcsNum
        % Extract IDs of cones connected to this RGC's RF center
        connectivityVector = full(squeeze(connectivityMatrix(:, RGCindex)));
        inputConeIDs = find(connectivityVector>0);
        rgcCenterConeInputsNum(RGCindex) = numel(inputConeIDs);
        if (rgcCenterConeInputsNum(RGCindex) == 0)
            error('RGC %d has zero inputs\n', RGCindex);
        end
        
        % Compute RGC RF center position computed by the summed inputs to the center
        rgcPositionsDegsFromConnectivity(RGCindex,:) = ...
            sum(bsxfun(@times, conePositionsDegs, connectivityVector),1)/sum(connectivityVector);
    end % RGCindex
    
    % Synthesize center/surround RF params both in retinal and visual domains
    synthesizedRFparams = RGCmodels.CronerKaplan.compute.synthesizedMidgetRGCRFparams(...
        rgcPositionsDegsFromConnectivity, rgcCenterConeInputsNum);

    % Indices for sparse matrices
    rgcIndicesVectorC = [];
    coneIndicesVectorC = [];
    weightsVectorC = [];
    rgcIndicesVectorS = [];
    coneIndicesVectorS = [];
    weightsVectorS = [];
    
    for RGCindex  = 1:rgcsNum
        fprintf('Generating weights for RGC %d of %d\n', RGCindex, rgcsNum);
        
        %
        % CENTER WEIGHTS
        %
        exclusiveConnectionsVector = full(squeeze(connectivityMatrix(:, RGCindex)));
        synthParams = struct(...
            'subregionPeakSensivity', synthesizedRFparams.retinal.centerPeakSensitivities(RGCindex) ...
            );
        [rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, centerWeightsForThisRGC] = computeSubregionConeWeights( ...
            'center', RGCindex, coneTypes, conePositionsDegs,  ...
            rgcPositionsDegsFromConnectivity(RGCindex,:), ... 
            exclusiveConnectionsVector, synthParams, ...
            rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC);
        
        %
        % SURROUND WEIGHTS
        %
        exclusiveConnectionsVector = [];
        synthParams = struct(...
            'subregionRadiusDegs', synthesizedRFparams.retinal.surroundCharacteristicRadiiDegs(RGCindex), ... 
            'subregionPeakSensivity', synthesizedRFparams.retinal.surroundPeakSensitivities(RGCindex) ...
            );
        [rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, surroundWeightsForThisRGC] = computeSubregionConeWeights( ...
            'surround', RGCindex, coneTypes, conePositionsDegs,  ...
            rgcPositionsDegsFromConnectivity(RGCindex,:), ...
            exclusiveConnectionsVector, synthParams, ...
            rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS);
        
        surroundToCenterIntegratedSensitivityRatioFromWeights = sum(surroundWeightsForThisRGC)/sum(centerWeightsForThisRGC);
        fprintf('Retinal surround/center integrated sensitivity ratio from weights: %2.2f\n', surroundToCenterIntegratedSensitivityRatioFromWeights);
    end % RGCindex
    
    % Form coneWeights struct containing sparse matrices with cone weights
    % to the center and the surround subregions
    sparseMatrixRows = conesNum;
    sparseMatrixCols = rgcsNum;
    coneWeights.surround = sparse(coneIndicesVectorS', rgcIndicesVectorS', weightsVectorS', sparseMatrixRows, sparseMatrixCols);
    coneWeights.center = sparse(coneIndicesVectorC', rgcIndicesVectorC', weightsVectorC', sparseMatrixRows, sparseMatrixCols);
end

function [rgcIndicesVector, coneIndicesVector, weightsVector, weights] = computeSubregionConeWeights( ...
            subregionName, RGCindex, coneTypes, conePositionsMicrons, rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVector, coneIndicesVector, weightsVector)
        
    if (strcmp(subregionName, 'center'))
        % Only accept cones with exclusive connections, and assign equal
        % weights to all of these inputs
        connectedConeIndices = find(exclusiveConnectionsVector>0);
        weights(1:numel(connectedConeIndices),1) = 1.0;
    else
        % This is a surround subregion, so we compute Gaussian cone weights
        weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, synthesizedParams.subregionRadiusDegs);
        connectedConeIndices = find(weights >= exp(-3));
        netWeightsIncludingScones = sum(weights(connectedConeIndices));
        % Find S-cones, and set their weight to 0.
        sConeIndices = find(coneTypes(connectedConeIndices) == mRGCmosaic.SCONE_ID);
        weights(connectedConeIndices(sConeIndices)) = 0.000001;
        netWeightsExcludingScones = sum(weights(connectedConeIndices));
        % Correct weights since S-cones do not connect
        sConeCorrection = netWeightsIncludingScones/netWeightsExcludingScones;
        weights = weights(connectedConeIndices) * sConeCorrection ;
    end
    
    % Multiply with subregion peak sensitivity
    weights = weights * synthesizedParams.subregionPeakSensivity;

    % Acummulate sparse matrix indices for the surround
    rgcIndicesVector = cat(1, rgcIndicesVector, repmat(RGCindex, [numel(connectedConeIndices) 1]));
    coneIndicesVector = cat(1, coneIndicesVector, connectedConeIndices);
    weightsVector = cat(1, weightsVector, weights);
end

function weights = gaussianConeWeights(conePositions, rgcPosition, rgcSubregionRadius)
    deltaPos = bsxfun(@minus, conePositions, rgcPosition);
    distanceFromCenterSquared = sum(deltaPos.^2,2);
    weights = exp(-distanceFromCenterSquared/rgcSubregionRadius^2);
end