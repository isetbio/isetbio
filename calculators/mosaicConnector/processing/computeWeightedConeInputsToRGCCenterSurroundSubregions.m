function [connectionMatrixCenter, connectionMatrixSurround, ...
    synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
        conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        midgetRGCconnectionMatrix, ...
        rgcMosaicPatchEccMicrons, rgcMosaicPatchSizeMicrons, ...
        deconvolutionOpticsParams, visualizePatchDeconvolutionModel, figExportsDir)
    
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Synthesize RF params using the Croner&Kaplan model
    synthesizedRFParams = synthesizeRFparams(conePositionsMicrons, coneSpacingsMicrons, midgetRGCconnectionMatrix, deconvolutionOpticsParams);
    
    % Visualize deconvolution model for this patch
    if (visualizePatchDeconvolutionModel)
        visualizeDeconvolutionModelForConnectedPatch(synthesizedRFParams, figExportsDir);
    end
    
    % Convert patch ecc and size from retinal microns to visual degrees
    synthesizedRFParams.patchEccDegs = WatsonRGCModel.rhoMMsToDegs(rgcMosaicPatchEccMicrons*1e-3);
    synthesizedRFParams.patchSizeDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(rgcMosaicPatchSizeMicrons, rgcMosaicPatchEccMicrons);
    
    % Retrieve the retinal synthesized params
    synthesizedRetinalParams = synthesizedRFParams.retinal;
    rgcCenterPositionMicrons = synthesizedRFParams.centerPositionMicrons;
    rgcEccDegs = synthesizedRFParams.eccDegs;
    
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
            'subregionRadiusMicrons', WatsonRGCModel.sizeDegsToSizeRetinalMicrons(synthesizedRetinalParams.centerRadiiDegs(iRGC), rgcEccDegs(iRGC)), ... 
            'subregionRadiusDegs',    synthesizedRetinalParams.centerRadiiDegs(iRGC), ...
            'subregionPeakSensivity', synthesizedRetinalParams.centerPeakSensitivities(iRGC) ...
            );
        [rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, centerWeightsForThisRGC, actualToModelAreaRatio] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons, coneSpacingsMicrons, rgcPositionMicrons, ... 
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, true, []);
        actualToModelAreaRatio
        %
        % SURROUND WEIGHTS
        %
        exclusiveConnectionsVector = [];
        synthesizedParams = struct(...
            'subregionRadiusMicrons', WatsonRGCModel.sizeDegsToSizeRetinalMicrons(synthesizedRetinalParams.surroundRadiiDegs(iRGC), rgcEccDegs(iRGC)), ... 
            'subregionRadiusDegs',    synthesizedRetinalParams.surroundRadiiDegs(iRGC), ...
            'subregionPeakSensivity', synthesizedRetinalParams.surroundPeakSensitivities(iRGC) ...
            );
        [rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, surroundWeightsForThisRGC] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons, coneSpacingsMicrons, rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, false, actualToModelAreaRatio);
        
       
        
       % Check synthesized and actual S/C integrated sensitivity ratio
        synthesizedSCintegratedSensitivityRatio = (synthesizedRetinalParams.surroundPeakSensitivities(iRGC)/synthesizedRetinalParams.centerPeakSensitivities(iRGC)) * ...
            (synthesizedRetinalParams.surroundRadiiDegs(iRGC)/synthesizedRetinalParams.centerRadiiDegs(iRGC))^2;
        actualSCintegratedSensitivityRatio = sum(surroundWeightsForThisRGC)/sum(centerWeightsForThisRGC);
        
        fprintf('Net weight ratio: %2.2f, model integrated sensitivity ratio: %2.2f, with %d center inputs\n', ...
            actualSCintegratedSensitivityRatio, ...
            synthesizedSCintegratedSensitivityRatio, ...
            numel(centerWeightsForThisRGC));
        
        
    end % iRGC
    
    % Form sparse matrices with weights
    connectionMatrixSurround = sparse(coneIndicesVectorS', rgcIndicesVectorS', weightsVectorS', conesNum, rgcsNum);
    connectionMatrixCenter = sparse(coneIndicesVectorC', rgcIndicesVectorC', weightsVectorC', conesNum, rgcsNum);

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


function [rgcIndicesVector, coneIndicesVector, weightsVector, weights, actualToModelAreaRatio] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons, coneSpacingsMicrons, rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVector, coneIndicesVector, weightsVector, printConnectionWeights, actualToModelAreaRatio)
        
    global SCONE_ID
    
    % Retrieve cone indices connected to the RGC surround (up to 3 sigma)
    weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, synthesizedParams.subregionRadiusMicrons);
    connectedConeIndices = find(weights >= exp(-4));
    
    if (~isempty(exclusiveConnectionsVector))
        % This is a center subregion
        coneIndicesWithExclusiveConnections = find(exclusiveConnectionsVector>0);
        connectedConeIndices = coneIndicesWithExclusiveConnections;
        
        inputsNum = numel(connectedConeIndices);
        modelRFdiameter1sigma = synthesizedParams.subregionRadiusMicrons*2;
        coneDiameters = coneSpacingsMicrons(connectedConeIndices)*0.7;
        meanConeDiameter = mean(coneDiameters);
        meanConeDiameter1sigma = meanConeDiameter/3;
        
        fprintf('%d-cone input center with model RF center diameter (1 sigma): %2.2f um, cone diameter: %2.2f um\n', ...
            inputsNum, modelRFdiameter1sigma, meanConeDiameter1sigma);

        % If center has less than or equal to 4 cones, do not apply
        % Gaussian cone weights. Instead make all inputs equal
        if (numel(coneIndicesWithExclusiveConnections)<=4)
            connectedConeIndices = coneIndicesWithExclusiveConnections;
            weights = weights * 0;
            weights(connectedConeIndices) = 1.0/numel(coneIndicesWithExclusiveConnections);
        end
    end
    
    % Find S-cones, and set their weight to 0.
    sConeIndices = find(coneTypes(connectedConeIndices) == SCONE_ID);
    weights(connectedConeIndices(sConeIndices)) = 0.000001;
    weights = weights(connectedConeIndices);
    
    % Find actual area occupied by the input cones
    [actualSensitivity, ~, actualSubregionAreaMicrons2At1Sigma] = ...
            computeActualSubregionIntegratedSensitivity(...
            conePositionsMicrons(connectedConeIndices,:), ...
            coneSpacingsMicrons(connectedConeIndices), ...
            coneTypes(connectedConeIndices), weights);
        
    
    modelSubregionSensitivity = pi * (synthesizedParams.subregionRadiusMicrons*3)^2;
    sensitivityCorrectionFactor = modelSubregionSensitivity / actualSensitivity;
    
    if (isempty(exclusiveConnectionsVector))
        % Boost surround
        %sensitivityCorrectionFactor  = sensitivityCorrectionFactor * actualToModelAreaRatio;
    else
        modelSubregionAreaMicrons2At1Sigma = pi*synthesizedParams.subregionRadiusMicrons^2;
        actualToModelAreaRatio = actualSubregionAreaMicrons2At1Sigma/modelSubregionAreaMicrons2At1Sigma;
    end
    
    
    % Multiply with subregion peak sensitivity
    weights = weights * sensitivityCorrectionFactor * synthesizedParams.subregionPeakSensivity;

    if (printConnectionWeights)
        fprintf('weights (%d, sum: %2.2f): ', numel(weights), sum(weights))
        for k = 1:numel(weights)
            fprintf('%2.0f ', weights(k))
        end
        fprintf('\n');
    end
    
    % Acummulate sparse matrix indices for the surround
    rgcIndicesVector = cat(1, rgcIndicesVector, repmat(iRGC, [numel(connectedConeIndices) 1]));
    coneIndicesVector = cat(1, coneIndicesVector, connectedConeIndices);
    weightsVector = cat(1, weightsVector, weights);
end


function [actualSensitivity, actualAreaDeg2, actualAreaMicrons2At1Sigma] = computeActualSubregionIntegratedSensitivity(...
          conePositionsMicrons, coneSpacingsMicrons, coneTypes, weights)
    global SCONE_ID
    actualSensitivity = 0;
    actualAreaDeg2 = 0;
    actualAreaMicrons2At1Sigma = 0;
    netWeight = 0;
    
    for k = 1:numel(coneSpacingsMicrons)
        if (coneTypes(k) ~= SCONE_ID)
            coneDiameterMicrons = coneSpacingsMicrons(k);
            coneApertureMicrons = 0.7*coneDiameterMicrons;
            
            coneEccMicrons = sqrt(sum(conePositionsMicrons(k,:).^2,2));
            coneApertureDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(coneApertureMicrons, coneEccMicrons);
            
            coneApertureFullRadiusDegs = 0.5*coneApertureDegs;
            coneApertureFullRadiusMicrons = 0.5*coneApertureMicrons;
            coneApertureSigmaRadiusMicrons = coneApertureFullRadiusMicrons/3;
            
            coneApertureAreaDegs2 = pi * coneApertureFullRadiusDegs^2;
            coneApertureAreaMicrons2 = pi * coneApertureFullRadiusMicrons^2;
            coneApertureAreaMicrons2At1Sigma = pi * coneApertureSigmaRadiusMicrons^2;
            
            actualSensitivity = actualSensitivity + weights(k) * coneApertureAreaMicrons2;
            netWeight = netWeight + weights(k);
            
            % For more than 1 cones, this area measurement may be
            % exagerrated as we are counting overlapping regions twice
            actualAreaDeg2 = actualAreaDeg2 + coneApertureAreaDegs2;
            actualAreaMicrons2At1Sigma = actualAreaMicrons2At1Sigma + coneApertureAreaMicrons2At1Sigma;
        end
    end
    
end

function weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSubregionRadiusMicrons)
    deltaPos = bsxfun(@minus, conePositionsMicrons, rgcPositionMicrons);
    weights = exp(-sum(deltaPos.^2,2)/rgcSubregionRadiusMicrons^2);
end