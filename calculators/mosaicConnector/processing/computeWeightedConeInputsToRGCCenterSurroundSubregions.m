function [connectionMatrixCenter, connectionMatrixSurround, ...
    synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
        conePositionsMicrons,  coneTypes, ...
        midgetRGCconnectionMatrix, ...
        rgcMosaicPatchEccMicrons, rgcMosaicPatchSizeMicrons, ...
        deconvolutionOpticsParams, visualizePatchDeconvolutionModel, figExportsDir)
    
    
    % Synthesize RF params using the Croner&Kaplan model
    synthesizedRFParams = synthesizeRFparams(conePositionsMicrons,  midgetRGCconnectionMatrix, deconvolutionOpticsParams);
    
    % Visualize deconvolution model for this patch
    if (visualizePatchDeconvolutionModel)
        visualizeDeconvolutionModelForConnectedPatch(synthesizedRFParams, figExportsDir);
    end
    
    % Convert patch ecc and size from retinal microns to visual degrees
    synthesizedRFParams.patchEccDegs = WatsonRGCModel.rhoMMsToDegs(rgcMosaicPatchEccMicrons*1e-3);
    synthesizedRFParams.patchSizeDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(rgcMosaicPatchSizeMicrons, rgcMosaicPatchEccMicrons);
    
    % Retrieve the retinal synthesized params
    synthesizedRetinalParams = synthesizedRFParams.retinal;
    rgcCenterPositionMicrons = synthesizedRFParams.rfCenterPositionMicrons;    
    
    % Compute weights of cone inputs to the RF center and to the RF
    % surround based on the midgetRGCconnectionMatrix and the
    % synthesizedRFparams
    
    conesNum = size(midgetRGCconnectionMatrix,1);
    rgcsNum = size(midgetRGCconnectionMatrix,2);
  
    % Indices for sparse matrices
    rgcIndicesVectorC = [];
    coneIndicesVectorC = [];
    weightsVectorC = [];
    rgcIndicesVectorS = [];
    coneIndicesVectorS = [];
    weightsVectorS = [];
    
    for iRGC = 1:rgcsNum
        fprintf('Generating weights for RGC %d of %d\n', iRGC, rgcsNum);
        
        % Get position of RGC in microns
        rgcPositionMicrons = rgcCenterPositionMicrons(iRGC,:);

        %
        % CENTER WEIGHTS
        %
        exclusiveConnectionsVector = full(squeeze(midgetRGCconnectionMatrix(:, iRGC)));
        synthesizedParams = struct(...
            'subregionRadiusMicrons', [], ... 
            'subregionPeakSensivity', synthesizedRetinalParams.centerPeakSensitivities(iRGC) ...
            );
        [rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, centerWeightsForThisRGC] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons,  rgcPositionMicrons, ... 
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC);
        %
        % SURROUND WEIGHTS
        %
        exclusiveConnectionsVector = [];
        synthesizedParams = struct(...
            'subregionRadiusMicrons', synthesizedRetinalParams.surroundCharacteristicRadiiMicrons(iRGC), ... 
            'subregionPeakSensivity', synthesizedRetinalParams.surroundPeakSensitivities(iRGC) ...
            );
        [rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, surroundWeightsForThisRGC] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons,  rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS);
        
       
        surroundToCenterIntegratedSensitivityRatio = sum(surroundWeightsForThisRGC)/sum(centerWeightsForThisRGC);
        
        fprintf('Retinal surround to center integrated sensitivity ratio: %2.2f\n', surroundToCenterIntegratedSensitivityRatio);
    end % iRGC
    
    % Form sparse matrices with weights
    sparseMatrixRows = conesNum;
    sparseMatrixCols = rgcsNum;
    %S(coneIndicesVectorS, rgcIndicesVectorS) = weightsVectorS
    connectionMatrixSurround = sparse(coneIndicesVectorS', rgcIndicesVectorS', weightsVectorS', sparseMatrixRows, sparseMatrixCols);
    connectionMatrixCenter = sparse(coneIndicesVectorC', rgcIndicesVectorC', weightsVectorC', sparseMatrixRows, sparseMatrixCols);

    % Find indices of RGCs within the [rgcMosaicPatchEccMicrons,
    % rgcMosaicPatchSizeMicrons] region, which have > 0 inputs to their
    % center
    rgcsNum = size(connectionMatrixCenter,2);
    rgcIndicesWithNonZeroInputs = [];
    for mRGCindex = 1:rgcsNum
        weights = full(connectionMatrixCenter(:,mRGCindex));
        centerIndices = find(weights>0);
        if (numel(centerIndices) > 0)
            rgcIndicesWithNonZeroInputs = cat(2,rgcIndicesWithNonZeroInputs, mRGCindex);
        end
    end % mRGCindex
    
    % Only keep weights for RGCs within the [rgcMosaicPatchEccMicrons, rgcMosaicPatchSizeMicrons] region
    connectionMatrixCenter = connectionMatrixCenter(:, rgcIndicesWithNonZeroInputs);
    connectionMatrixSurround = connectionMatrixSurround(:, rgcIndicesWithNonZeroInputs);
end


function [rgcIndicesVector, coneIndicesVector, weightsVector, weights] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons, rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVector, coneIndicesVector, weightsVector)
        
    global SCONE_ID
    
    if (~isempty(exclusiveConnectionsVector))
        % This is a center subregion
        coneIndicesWithExclusiveConnections = find(exclusiveConnectionsVector>0);
        % Only accept cones with exclusive connections
        connectedConeIndices = coneIndicesWithExclusiveConnections;
        weights(1:numel(connectedConeIndices),1) = 1.0;
    else
        % This is a surround subregion
        weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, synthesizedParams.subregionRadiusMicrons);
        connectedConeIndices = find(weights >= exp(-4));
        % Find S-cones, and set their weight to 0.
        sConeIndices = find(coneTypes(connectedConeIndices) == SCONE_ID);
        weights(connectedConeIndices(sConeIndices)) = 0.000001;
        weights = weights(connectedConeIndices);
    end
    
    % Multiply with subregion peak sensitivity
    weights = weights * synthesizedParams.subregionPeakSensivity;

    % Acummulate sparse matrix indices for the surround
    rgcIndicesVector = cat(1, rgcIndicesVector, repmat(iRGC, [numel(connectedConeIndices) 1]));
    coneIndicesVector = cat(1, coneIndicesVector, connectedConeIndices);
    weightsVector = cat(1, weightsVector, weights);
end

function weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSubregionRadiusMicrons)
    deltaPos = bsxfun(@minus, conePositionsMicrons, rgcPositionMicrons);
    weights = exp(-sum(deltaPos.^2,2)/rgcSubregionRadiusMicrons^2);
end