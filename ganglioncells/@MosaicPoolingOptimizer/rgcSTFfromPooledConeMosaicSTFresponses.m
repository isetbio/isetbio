function theSTFdata = rgcSTFfromPooledConeMosaicSTFresponses(obj, pooledConeIndicesAndWeights, visualRcDegs)

    % Compute the pooled center signals for all sfs, all orientations, all spatial phases
    conesNum = numel(pooledConeIndicesAndWeights.centerConeIndices);
    centerConeWeights = reshape(pooledConeIndicesAndWeights.centerConeWeights(:), [1 1 1 conesNum]);
    centerConeResponses = obj.inputConeMosaicVisualSTFdata.responseModulations(:,:,:, pooledConeIndicesAndWeights.centerConeIndices);
    netCenterResponses = sum(bsxfun(@times, centerConeResponses, centerConeWeights), 4);

    % Compute the pooled surround signals
    conesNum = numel(pooledConeIndicesAndWeights.surroundConeWeights);
    surroundConeWeights = reshape(pooledConeIndicesAndWeights.surroundConeWeights(:), [1 1 1 conesNum]);
    surroundConeResponses = obj.inputConeMosaicVisualSTFdata.responseModulations(:,:,:, pooledConeIndicesAndWeights.surroundConeIndices);
    netSurroundResponses = sum(bsxfun(@times, surroundConeResponses, surroundConeWeights), 4);

    % Composite responses (center-surround) of the current modelRGC (defined by the current pooledConeIndicesAndWeights)
    theModelRGCResponsesForAllOrientationsAndSpatialFrequencies = netCenterResponses - netSurroundResponses;

    % The STF of the current model RGC
    orientationsTested = obj.inputConeMosaicVisualSTFdata.orientationsTested;
    spatialFrequenciesTested = obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested;
    theModelRGCoptimalSTF = MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
                    orientationsTested, spatialFrequenciesTested, ...
                    theModelRGCResponsesForAllOrientationsAndSpatialFrequencies);


    % Fit the visual STF with a DoG model
    [theFittedDoGModelParams, theFittedDoGModelToTheVisualSTF] = MosaicPoolingOptimizer.fitDifferenceOfGaussiansToSTF(...
                      obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested, ...
                      theModelRGCoptimalSTF, ...
                      visualRcDegs, ...
                      visualRcDegs*[0.5 1 1.5], ...
                      obj.multiStartsNumDoGFit);

    theSTFdata = struct;
    theSTFdata.spatialFrequencySupport = obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested;
    theSTFdata.visualSTF = theModelRGCoptimalSTF;
    theSTFdata.fittedDoGModelParams = theFittedDoGModelParams;
    theSTFdata.fittedDoGModelToVisualSTF = theFittedDoGModelToTheVisualSTF;
    theSTFdata.fittedRcDegs = theFittedDoGModelParams.finalValues(4);
    theSTFdata.fittedDoGModelRsRcRatio = theFittedDoGModelParams.finalValues(3);
    theSTFdata.fittedDoGModelSCIntSensRatio = theFittedDoGModelParams.finalValues(2) * (theFittedDoGModelParams.finalValues(3))^2;
    
end


