function [connectionMatrixCenter, connectionMatrixSurround, ...
    synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
        conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
        midgetRGCconnectionMatrix, ...
        rgcMosaicPatchEccMicrons, rgcMosaicPatchSizeMicrons, ...
        deconvolutionOpticsParams)
    
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Synthesize RF params using the Croner&Kaplan model
    synthesizedRFParams = synthesizeRFparams(conePositionsMicrons, coneSpacingsMicrons, midgetRGCconnectionMatrix, deconvolutionOpticsParams);

    % Convert patch ecc and size from retinal microns to visual degrees
    synthesizedRFParams.patchEccDegs = WatsonRGCModel.rhoMMsToDegs(rgcMosaicPatchEccMicrons*1e-3);
    synthesizedRFParams.patchSizeDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(rgcMosaicPatchSizeMicrons, rgcMosaicPatchEccMicrons);
    
    rgcParams = synthesizedRFParams.retinal;
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
        
        % Get index of RGC in full mosaic
        rgcIndex = iRGC;
        
        % Get position of RGC in microns
        rgcPositionMicrons = rgcCenterPositionMicrons(iRGC,:);
        
        % Get radius of RGC center in microns
        rgcCenterRadiusDegs = rgcParams.centerRadiiDegs(iRGC);
        rgcCenterRadiusMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(rgcCenterRadiusDegs, rgcEccDegs(iRGC));
        
        % Get peak sensitivity of RGC center
        rgcCenterPeakSensivity = rgcParams.centerPeakSensitivities(iRGC);
        
        % Get radius of RGC surround in microns
        rgcSurroundRadiusDegs = rgcParams.surroundRadiiDegs(iRGC);
        rgcSurroundRadiusMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(rgcSurroundRadiusDegs, rgcEccDegs(iRGC));
        
        % Get peak sensitivity of RGC surround
        rgcSurroundPeakSensivity = rgcParams.surroundPeakSensitivities(iRGC);
        
        % Retrieve cone indices connected to the RGC surround, 99%
        weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSurroundRadiusMicrons);
        coneIndicesConnectedToSurround = find(weights >= 0.01);

        % Find S-cones, and set their weight to 0, as H1 horizontal cells
        % (which make the surrounds of mRGCs) do not contact S-cones.
        sConeIndices = find(coneTypes(coneIndicesConnectedToSurround) == SCONE_ID);
        weights(coneIndicesConnectedToSurround(sConeIndices)) = 0.000001;
        weightsS = weights(coneIndicesConnectedToSurround);
        weightsS = weightsS * rgcSurroundPeakSensivity;
        
        % Acummulate sparse matrix indices for the surround
        rgcIndicesVectorS = cat(1, rgcIndicesVectorS, repmat(rgcIndex, [numel(coneIndicesConnectedToSurround) 1]));
        coneIndicesVectorS = cat(1, coneIndicesVectorS, coneIndicesConnectedToSurround);
        weightsVectorS = cat(1, weightsVectorS, weightsS);
        
        % Retrieve cone indices connected to the RGC center. These have
        % been determined by the local topology of the cone mosaic and the
        % mRGCRF mosaic, in a way that avoid all overlap. That is each cone
        % is connected to only one RGC. 
        connectivityVector = full(squeeze(midgetRGCconnectionMatrix(:, rgcIndex)));
        coneIndicesConnectedExclusivelyToCenter = find(connectivityVector>0);
        
        % We can use the RGC center radius (from the Croner&Kaplan model)
        % to pool more cones into the center (but these will be providing
        % most of their output to other RGCs). These additional cones will
        % result in some RF center overlap
        acceptSharedConeInputWithinRadiusMatchingTheCronerKaplanModel = true;
        if (acceptSharedConeInputWithinRadiusMatchingTheCronerKaplanModel)
            % Retrieve cone indices connected to the RGC center, 90%
            weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcCenterRadiusMicrons);
            coneIndicesNearRFcenter = find((weights >= 0.1)&((coneTypes==LCONE_ID)|(coneTypes==MCONE_ID)));
        else
            coneIndicesNearRFcenter = [];
        end
        
        allConeIndicesConnectedToCenter = unique([coneIndicesConnectedExclusivelyToCenter(:); coneIndicesNearRFcenter(:)]);
        coneIndicesConnectedToCenterProvidingSharedInput = setdiff(allConeIndicesConnectedToCenter, coneIndicesConnectedExclusivelyToCenter);
        
        if (~isempty(coneIndicesConnectedToCenterProvidingSharedInput))
            fprintf(' RGC %d center: %d cones provide exclusive input and %d cones provide shared input\n', ...
                iRGC, numel(coneIndicesConnectedExclusivelyToCenter), numel(coneIndicesConnectedToCenterProvidingSharedInput));
        end
        
        % Gaussian weights for all cones
        weights = gaussianConeWeights(conePositionsMicrons(allConeIndicesConnectedToCenter,:), rgcPositionMicrons, rgcCenterRadiusMicrons);
        
        % Except for the origiMaximum weights (1) for cones 
        for iConeExclusive = 1:numel(coneIndicesConnectedExclusivelyToCenter)
            idx = find(allConeIndicesConnectedToCenter == coneIndicesConnectedExclusivelyToCenter(iConeExclusive));
            weights(idx) = 1.0;
        end
        
        % Multiple by center's peak sensitivity
        weightsC = weights * rgcCenterPeakSensivity;
        
        % Adjust center weights to ensure the achieved integrated sensitivity ratio matches the desired one
        desiredSurroundToCenterIntegratedSensitivityRatio = rgcSurroundPeakSensivity/rgcCenterPeakSensivity * (rgcSurroundRadiusMicrons/rgcCenterRadiusMicrons)^2;
        actualSurroundToCenterIntegratedSensitivityRatio = sum(weightsS)/sum(weightsC);

        sensitivityCorrectionFactor = desiredSurroundToCenterIntegratedSensitivityRatio/actualSurroundToCenterIntegratedSensitivityRatio;   
        weightsC = weightsC / sensitivityCorrectionFactor;
        
        % Acummulate sparse matrix indices for the surround
        rgcIndicesVectorC = cat(1, rgcIndicesVectorC, repmat(rgcIndex, [numel(allConeIndicesConnectedToCenter) 1]));
        coneIndicesVectorC = cat(1, coneIndicesVectorC, allConeIndicesConnectedToCenter);
        weightsVectorC = cat(1, weightsVectorC, weightsC);
    end
    
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


function weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSubregionRadiusMicrons)
    deltaPos = bsxfun(@minus, conePositionsMicrons, rgcPositionMicrons);
    weights = exp(-sum(deltaPos.^2,2)/rgcSubregionRadiusMicrons^2);
end