function compute(obj, gridNodeIndex, whichConeType, coneMosaicSTFresponsesFileName, optimizedRGCpoolingObjectsFileName, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('targetSurroundToCenterRcRatio', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('targetSurroundToCenterIntegratedSensitivityRatio', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('displayFittingProgress', false, @islogical);
    p.addParameter('multiStartsNumDoGFit', 64, @isscalar);
    p.addParameter('multiStartsNumRetinalPooling', 8, @isscalar);
    p.addParameter('rmseWeightForRsRcResidual', 1.0, @isscalar);
    p.addParameter('rmseWeightForSCintSensResidual', 1.0, @isscalar);
    p.addParameter('retinalRFmodelParams', MosaicPoolingOptimizer.defaultRetinalRFmodelParams, @isstruct);

    p.parse(varargin{:});

    obj.multiStartsNumDoGFit = p.Results.multiStartsNumDoGFit;
    obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
    obj.rmseWeightForRsRcResidual = abs(p.Results.rmseWeightForRsRcResidual);
    obj.rmseWeightForSCintSensResidual = abs(p.Results.rmseWeightForSCintSensResidual);
    obj.retinalRFmodelParams = p.Results.retinalRFmodelParams;

    displayFittingProgress = p.Results.displayFittingProgress;
    
    % Default targetSurroundToCenterRcRatio
    targetSurroundToCenterRcRatio = obj.visualSTFSurroundToCenterRcRatioGrid(gridNodeIndex);
    if (~isempty(p.Results.targetSurroundToCenterRcRatio))
        fprintf(2, 'Overriding the default target Rs/Rc ratio of %2.2f with %2.2f\n', ...
             targetSurroundToCenterRcRatio, p.Results.targetSurroundToCenterRcRatio);
        targetSurroundToCenterRcRatio = p.Results.targetSurroundToCenterRcRatio;
    end

    % Default targetSurroundToCenterIntegratedSensitivityRatio
    targetSurroundToCenterIntegratedSensitivityRatio = obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid(gridNodeIndex);
    if (~isempty(p.Results.targetSurroundToCenterIntegratedSensitivityRatio))
        fprintf(2, 'Overriding the default Rs/Rc ratio (%2.2f) with %2.2f\n', ...
             targetSurroundToCenterIntegratedSensitivityRatio, p.Results.targetSurroundToCenterIntegratedSensitivityRatio);
        targetSurroundToCenterIntegratedSensitivityRatio = p.Results.targetSurroundToCenterIntegratedSensitivityRatio;
    end

    targetVisualSTFparams = struct(...
        'surroundToCenterRcRatio', targetSurroundToCenterRcRatio, ...
        'surroundToCenterIntegratedSensitivityRatio', targetSurroundToCenterIntegratedSensitivityRatio ...
        );

    % Optimized RGCpooling object filename
    optimizedRGCpoolingObjectsFileNameForThisNode = ...
        strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));
    fprintf('Will save optimized object at %s\n', optimizedRGCpoolingObjectsFileNameForThisNode)

    % Load the precomputed cone mosaic STF responses
    obj.loadConeMosaicVisualSTFresponses(coneMosaicSTFresponsesFileName);

    % Deal with initialRetinalConePoolingParamds
    if (isfile(optimizedRGCpoolingObjectsFileNameForThisNode))
        fprintf('<<<<< Loading previously computed model\n');
        load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
            'theLconeRFcomputeStruct', 'theMconeRFcomputeStruct');
        initialRetinalLconePoolingParams = theLconeRFcomputeStruct.retinalConePoolingParams.finalValues;
        initialRetinalMconePoolingParams = theMconeRFcomputeStruct.retinalConePoolingParams.finalValues;
    else
        fprintf('Did not find a previously computed model. Starting with default params.\n')
        initialRetinalLconePoolingParams = [];
        initialRetinalMconePoolingParams = [];
    end

    % Optimize the L-center RGC RF pooling
    if (ismember(cMosaic.LCONE_ID, whichConeType))
        % Compute theLconeRFcomputeStruct
        fprintf('>>>>> Recomputing the LconeComputeStruct\n');
        
        LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
        figNo = 1000 + gridNodeIndex;
        figTitle = sprintf('grid no %d of %d L-cone center RGC %d', gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), LconeRGCindex);
        theLconeRFcomputeStruct = obj.optimizeSurroundConePooling(...
            LconeRGCindex, targetVisualSTFparams, ...
            initialRetinalLconePoolingParams, displayFittingProgress, figNo, figTitle);
    end

    % Optimize the M-center RGC RF pooling
    if (ismember(cMosaic.MCONE_ID, whichConeType))
        fprintf('>>>> Recomputing the MconeComputeStruct\n');
        MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
        figNo = 2000 + gridNodeIndex;
        figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
            gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);
    
        theMconeRFcomputeStruct = obj.optimizeSurroundConePooling(...
            MconeRGCindex, targetVisualSTFparams, ...
            initialRetinalMconePoolingParams, displayFittingProgress, figNo, figTitle);
    end

    % Saved computed object
    fprintf('\n\nSaved computeStructs for grid node %d to %s\n\n',gridNodeIndex, optimizedRGCpoolingObjectsFileNameForThisNode);
    save(optimizedRGCpoolingObjectsFileNameForThisNode, ...
        'theLconeRFcomputeStruct', ...
        'theMconeRFcomputeStruct');
end



