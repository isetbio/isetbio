function compute(obj, gridNodeIndex, coneMosaicSTFresponsesFileName, optimizedRGCpoolingObjectsFileName, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('targetSurroundToCenterRcRatio', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('targetSurroundToCenterIntegratedSensitivityRatio', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('displayFittingProgress', false, @islogical);
    p.addParameter('multiStartsNumDoGFit', 64, @isscalar);
    p.addParameter('multiStartsNumRetinalPooling', 8, @isscalar);

    p.parse(varargin{:});

    obj.multiStartsNumDoGFit = p.Results.multiStartsNumDoGFit;
    obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
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

    
    % Optimize the L-center RGC RF pooling
    LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
    figNo = 1000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d L-cone center RGC %d', gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), LconeRGCindex);
    theLconeRFcomputeStruct = obj.optimizeSurroundConePooling(...
        LconeRGCindex, targetVisualSTFparams, ...
        displayFittingProgress, figNo, figTitle);

    % Fit the M-center RGC RF pooling
    MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
    figNo = 2000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);

    theMconeRFcomputeStruct = obj.optimizeSurroundConePooling(...
        MconeRGCindex, targetVisualSTFparams, ...
        displayFittingProgress, figNo, figTitle);


    % Saved computed object
    fprintf('\n\nSaved computeStructs for grid node %d to %s\n\n',gridNodeIndex, optimizedRGCpoolingObjectsFileNameForThisNode);
    save(optimizedRGCpoolingObjectsFileNameForThisNode, ...
        'theLconeRFcomputeStruct', ...
        'theMconeRFcomputeStruct');
end




