function inspect(obj, gridNodeIndex, optimizedRGCpoolingObjectsFileName, varargin)
    % Parse optional input
%     p = inputParser;
%     p.addParameter('visualize', [], @(x)(isempty(x)||(isscalar(x))));
%     p.parse(varargin{:});

    % Optimized RGCpooling object filename
    optimizedRGCpoolingObjectsFileNameForThisNode = ...
        strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

    fprintf('Loading optimized L- and M-cone RF compute structs from %s\n', optimizedRGCpoolingObjectsFileNameForThisNode);
    load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
        'theLconeRFcomputeStruct', ...
        'theMconeRFcomputeStruct');

    LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
    figNo = 10000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d L-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), LconeRGCindex);
    hFig = inspectConeSpecificRFcomputeStruct(figNo, figTitle,theLconeRFcomputeStruct);

    MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
    figNo = 20000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);
    hFig = inspectConeSpecificRFcomputeStruct(figNo, figTitle,theMconeRFcomputeStruct);

end

function hFig = inspectConeSpecificRFcomputeStruct(figNo, figTitle, ...
    theConeSpecificRFcomputeStruct)

    % Retrieve the saved data
    targetVisualSTFparams = theConeSpecificRFcomputeStruct.theTargetSTFparams;
    theFinalSTFdata = theConeSpecificRFcomputeStruct.theAchievedSTFdata;
    retinalConePoolingParams = theConeSpecificRFcomputeStruct.retinalConePoolingParams;
    theFinalPooledConeIndicesAndWeights = theConeSpecificRFcomputeStruct.theFinalPooledConeIndicesAndWeights;
    rmseSequence = [];

    hFig = MosaicPoolingOptimizer.visualizeOptimizationProgress(...
                figNo, figTitle, ...
                targetVisualSTFparams, ...
                theFinalSTFdata, retinalConePoolingParams, ...
                theFinalPooledConeIndicesAndWeights, ...
                rmseSequence);

end

