function run()

    % Ask user which mosaic to use
    mosaicEcc = MosaicPoolingOptimizer.chooseMosaicToUse();

    % Get mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEcc);

    % Get operation to perform
    operationSetToPerformContains = MosaicPoolingOptimizer.operationsMenu(mosaicParams);

    % Perform the generateRGCMosaic operation
    if (operationSetToPerformContains.generateCenterConnectedRGCMosaic)
        if (isfield(mosaicParams, 'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs'))
            MosaicPoolingOptimizer.performGenerateCenterConnectedRGCMosaicOp(mosaicParams, ...
                'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', mosaicParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs);
        else
            MosaicPoolingOptimizer.performGenerateCenterConnectedRGCMosaicOp(mosaicParams);
        end
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
        % This controls the visualized spatial support range
        visualizedSpatialSupportRangeArcMin = 24;

        % This controls the visualize spatial frequency range.
        % [] for full range, or enter a specific range, i.e. [0.1 30]
        visualizedSpatialFrequencyRange = [0.1 30]; %[0.1 30];

        MosaicPoolingOptimizer.performInspectOptimizedSurroundConePoolingModelsOp(...
            mosaicParams, ...
            'tickSeparationArcMin', ceil(visualizedSpatialSupportRangeArcMin/6), ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);

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
        mosaicEccsToInclude = [0.0 2.5 7.0 -10.0];
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsToInclude);
        return;
    end

    % Perform the computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams, ...
            'maxSFcyclesPerDegree', 100, ...
            'stimSizeDegs', 0.3, ....
            'posIncrementDegs', max(mosaicParams.sizeDegs), ...
            'reComputeInputConeMosaicSubspaceRFmappingResponses', true, ...
            'reComputeMRGCMosaicSubspaceRFmappingResponses', false, ...
            'reComputeRFs', false ...
            );
        return;
    end

    % Perform the visualizeVisualRFmapForTargetRGC operation
    if (operationSetToPerformContains.visualizeVisualRFmapForTargetRGC)
        MosaicPoolingOptimizer.performVisualizeVisualRFmapForTargetRGC(mosaicParams);
        return;
    end

end


