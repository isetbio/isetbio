function animateModelConvergence(obj, gridNodeIndex, optimizedRGCpoolingObjectsFileName, varargin)
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

    theLconeRFcomputeStruct
    LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
    figNo = 1;
    figTitle = '';
    pdfFilename = '';
    animateConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
        obj.theRGCMosaic, ...
        theLconeRFcomputeStruct);


end


function animateConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
    theRGCMosaic, theConeSpecificRFcomputeStruct)

    % Retrieve the saved data
    theConeSpecificRFcomputeStruct
    targetVisualSTFparams = theConeSpecificRFcomputeStruct.theTargetSTFparams;
    theFinalSTFdata = theConeSpecificRFcomputeStruct.theAchievedSTFdata;
    retinalConePoolingParams = theConeSpecificRFcomputeStruct.retinalConePoolingParams;
    retinalConePoolingModel = theConeSpecificRFcomputeStruct.modelConstants.retinalConePoolingModel;
    theFinalPooledConeIndicesAndWeights = theConeSpecificRFcomputeStruct.theFinalPooledConeIndicesAndWeights;

    theConeSpecificRFcomputeStruct.theTargetSTFparams
    modelConstants = theConeSpecificRFcomputeStruct.modelConstants;
    theConeSpecificRFcomputeStruct.retinalConePoolingParams
    
    visualizeComponents = false;
    [modelConstants, retinalConePoolingParams, visualRcDegs] = ...
        obj.computeOptimizationComponents(theRGCindex, visualizeComponents);

    for iteration = 1:size(theConeSpecificRFcomputeStruct.rmseSequence,1)

        pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(...
            modelConstants, currentRetinalPoolingParamValues);

        theCurrentSTFdata = obj.rgcSTFfromPooledConeMosaicSTFresponses(...
            pooledConeIndicesAndWeights, visualRcDegs);

        rmse = theConeSpecificRFcomputeStruct.rmseSequence(iteration,:);
        RsRcRatioResidual = rmse(2);
        SCintSensRatioResidual = rmse(3);

        theConeSpecificRFcomputeStruct.theTargetSTFparams.theConeSpecificRFcomputeStruct.theTargetSTFparams
    end

end

