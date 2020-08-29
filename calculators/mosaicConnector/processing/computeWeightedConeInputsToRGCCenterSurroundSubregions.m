function [connectionMatrixCenter, connectionMatrixSurround, ...
    synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
        conePositionsMicrons,  coneTypes, ...
        midgetRGCconnectionMatrix, ...
        rgcMosaicPatchEccMicrons, rgcMosaicPatchSizeMicrons, ...
        deconvolutionOpticsParams, visualizeSynthesizedParams, figExportsDir)
    
    
    % Synthesize RF params using the Croner&Kaplan model
    synthesizedRFParams = synthesizeRFparams(conePositionsMicrons,  midgetRGCconnectionMatrix, deconvolutionOpticsParams);
    
    % Visualize deconvolution model for this patch
    if (visualizeSynthesizedParams)
        plotlabOBJ = setupPlotLab(0, 18, 10);
        visualizeSynthesizedCenterSurroundProperties(300, synthesizedRFParams, 'visual', plotlabOBJ, 'SynthesizedVisualParams', figExportsDir);
        visualizeSynthesizedCenterSurroundProperties(301, synthesizedRFParams, 'retinal', plotlabOBJ, 'SynthesizedRetinalParams', figExportsDir);
        setupPlotLab(-1);
    end
    
    % Convert patch ecc and size from retinal microns to visual degrees
    synthesizedRFParams.patchEccDegs = WatsonRGCModel.rhoMMsToDegs(rgcMosaicPatchEccMicrons*1e-3);
    synthesizedRFParams.patchSizeDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(rgcMosaicPatchSizeMicrons, rgcMosaicPatchEccMicrons);
    
    % Retrieve the retinal synthesized params
    synthesizedRetinalParams = synthesizedRFParams.retinal;
    rgcCenterPositionMicrons = synthesizedRFParams.rfCenterPositionMicrons;    
    
    % Compute weights of cone inputs to the RF center and to the RF
    % surround based on the midgetRGCconnectionMatrix and the synthesizedRFparams
    
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
        [rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, centerWeightsForThisRGC, overallGain] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons,  rgcPositionMicrons, ... 
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorC, coneIndicesVectorC, weightsVectorC, []);
        %
        % SURROUND WEIGHTS
        %
        exclusiveConnectionsVector = [];
        synthesizedParams = struct(...
            'subregionRadiusMicrons', synthesizedRetinalParams.surroundCharacteristicRadiiMicrons(iRGC), ... 
            'subregionPeakSensivity', synthesizedRetinalParams.surroundPeakSensitivities(iRGC) ...
            );
        [rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, surroundWeightsForThisRGC, ~] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons,  rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVectorS, coneIndicesVectorS, weightsVectorS, overallGain);
        
        surroundToCenterIntegratedSensitivityRatioFromWeights = sum(surroundWeightsForThisRGC)/sum(centerWeightsForThisRGC);
        fprintf('Retinal surround/center integrated sensitivity ratio from weights: %2.2f\n', surroundToCenterIntegratedSensitivityRatioFromWeights);
    end % iRGC
    
    % Form sparse matrices with weights
    sparseMatrixRows = conesNum;
    sparseMatrixCols = rgcsNum;
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


function [rgcIndicesVector, coneIndicesVector, weightsVector, weights, overallGain] = computeWeightsToSubregion( ...
            iRGC, coneTypes, conePositionsMicrons, rgcPositionMicrons, ...
            exclusiveConnectionsVector, synthesizedParams, ...
            rgcIndicesVector, coneIndicesVector, weightsVector, overallGain)
        
    global SCONE_ID
    
    if (~isempty(exclusiveConnectionsVector))
        % This is a center subregion
        coneIndicesWithExclusiveConnections = find(exclusiveConnectionsVector>0);
        % Only accept cones with exclusive connections
        connectedConeIndices = coneIndicesWithExclusiveConnections;
        weights(1:numel(connectedConeIndices),1) = 1.0;
        if (isempty(overallGain))
            overallGain = 1.0/numel(connectedConeIndices);
        end
        overallGain = 1;
    else
        % This is a surround subregion
        weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, synthesizedParams.subregionRadiusMicrons);
        connectedConeIndices = find(weights >= exp(-5));
        netWeightsIncludingScones = sum(weights(connectedConeIndices));
        % Find S-cones, and set their weight to 0.
        sConeIndices = find(coneTypes(connectedConeIndices) == SCONE_ID);
        weights(connectedConeIndices(sConeIndices)) = 0.000001;
        netWeightsExcludingScones = sum(weights(connectedConeIndices));
        % Correct weights since S-cones do not connect
        sConeCorrection = netWeightsIncludingScones/netWeightsExcludingScones;
        weights = weights(connectedConeIndices) * sConeCorrection ;
    end
    
    % Multiply with subregion peak sensitivity
    weights = weights * synthesizedParams.subregionPeakSensivity * overallGain;

    % Acummulate sparse matrix indices for the surround
    rgcIndicesVector = cat(1, rgcIndicesVector, repmat(iRGC, [numel(connectedConeIndices) 1]));
    coneIndicesVector = cat(1, coneIndicesVector, connectedConeIndices);
    weightsVector = cat(1, weightsVector, weights);
end

function weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSubregionRadiusMicrons)
    deltaPos = bsxfun(@minus, conePositionsMicrons, rgcPositionMicrons);
    distanceFromCenterSquared = sum(deltaPos.^2,2);
    weights = exp(-distanceFromCenterSquared/rgcSubregionRadiusMicrons^2);
end

function visualizeSynthesizedCenterSurroundProperties(figNo, synthesizedRFParams, dataSet, plotlabOBJ, pdfFileName, figExportsDir)
    
    % Get data
    switch (dataSet)
        case 'visual'
            rfParams = synthesizedRFParams.visual;
      
        case 'retinal'
            rfParams = synthesizedRFParams.retinal;
        otherwise
            error('Can only visualize params of the visual dataset. dataSet: ''%s''.', dataSet)
    end
    
    rfEccRadiusDegs = synthesizedRFParams.rfEccRadiusDegs;
    visualizeRFparamsForConnectedPatch(figNo, sprintf('%s sparams', dataSet), ...
        [], ...
        rfEccRadiusDegs, ...
        rfParams.centerCharacteristicRadiiDegs, rfParams.surroundCharacteristicRadiiDegs, ...
        rfParams.centerPeakSensitivities, rfParams.surroundPeakSensitivities, ...
        sprintf('%s__%s',pdfFileName, dataSet), ...
        figExportsDir, plotlabOBJ);
    
end

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 