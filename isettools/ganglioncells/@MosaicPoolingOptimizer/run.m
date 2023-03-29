function run()
    
    % Eccentricities of optimized mosaics
    mosaicEcc = 2.5;
    mosaicEcc = 7.0;
    mosaicEcc = -10.0;

    % Get mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEcc);

    % Get operation to perform
    operationSetToPerformContains = MosaicPoolingOptimizer.operationsMenu(mosaicParams);

    % Perform the generateRGCMosaic operation
    if (operationSetToPerformContains.generateCenterConnectedRGCMosaic)
        MosaicPoolingOptimizer.performGenerateCenterConnectedRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the visualizeRGCMosaic operation
    if (operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs)
        MosaicPoolingOptimizer.performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams);
        return;
    end

    % Perform the visualizePSFsWithinRGCMosaic operation
    if (operationSetToPerformContains.visualizePSFsWithinRGCMosaic)
        MosaicPoolingOptimizer.performVisualizePSFsWithinRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the computeInputConeMosaicSTFresponses operation
    if (operationSetToPerformContains.computeInputConeMosaicSTFresponses)
        MosaicPoolingOptimizer.performComputeInputConeMosaicSTFresponsesOp(mosaicParams);
        return;
    end

    % Perfom the optimizeSurroundConePoolingModels operation
    if (operationSetToPerformContains.optimizeSurroundConePoolingModels)
        
        % Fitting options
        multiStartsNumRetinalPooling = 20; % 12
        multiStartsNumDoGFit = 128;

        % More weight for matching the Rs/Rc ratio
        rmseWeightForRsRcResidual = 2.0;
        rmseWeightForSCintSensResidual = 1.0;

        MosaicPoolingOptimizer.performOptimizeSurroundConePoolingModelsOp(mosaicParams, ...
            'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
            'multiStartsNumDoGFit', multiStartsNumDoGFit, ...
            'rmseWeightForRsRcResidual', rmseWeightForRsRcResidual , ...
            'rmseWeightForSCintSensResidual', rmseWeightForSCintSensResidual);

        return;
    end


    % Perfom the inspectOptimizedSurroundConePoolingModels operation
    if (operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels)
        MosaicPoolingOptimizer.performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams);
        return;
    end

    % Perform the generateComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performGenerateComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performComputeVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performFitVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the visualizeConePoolingRFmapAndVisualSTFforTargetRGC operation
    if (operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC)
        MosaicPoolingOptimizer.performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics)
        mosaicEccsToInclude = [2.5 7.0];
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsToInclude);
        return;
    end

    % Perform the computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeVisualRFmapForTargetRGC operation
    if (operationSetToPerformContains.visualizeVisualRFmapForTargetRGC)
        MosaicPoolingOptimizer.performVisualizeVisualRFmapForTargetRGC(mosaicParams);
        return;
    end

end


