function run()

    % Ask user which mosaic to use
    mosaicEcc = MosaicPoolingOptimizer.chooseMosaicToUse();

    % Get mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEcc);

    % This controls the visualized spatial support range
    tickSeparationArcMin = 4;
    normalizedPeakSurroundSensitivity = 0.4;
    visualizedSpatialFrequencyRange = [0.1 100];

    %MosaicPoolingOptimizer.performComputeVisuallyProjectedAnatomicalConeRcsByFittingRFmap([]);

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
        
        identifyPooledCones = true;
        identifyInputCones = true;
        plotRFoutlines = true;
        labelRetinalMeridians  = true;
        backgroundColor = [1 1 1];


        MosaicPoolingOptimizer.performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones', identifyInputCones, ...
            'plotRFoutlines', plotRFoutlines, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'backgroundColor', backgroundColor);
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

    % Perform the computeInputConeMosaicVisuallyProjectedCharacteristicRadii operation
    if (operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii)
        MosaicPoolingOptimizer.performComputeInputConeMosaicVisuallyProjectedRcDegsOp(mosaicParams);
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

        MosaicPoolingOptimizer.performInspectOptimizedSurroundConePoolingModelsOp(...
            mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
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

    % Perform the adjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata
    if (operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs)
        MosaicPoolingOptimizer.performAdjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata(mosaicParams)
        return;
    end


    % Perform the visualizeConePoolingRFmapAndVisualSTFforTargetRGC operation
    if (operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC)
        
        MosaicPoolingOptimizer.performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics)
        mosaicEccsToInclude = [0.0 2.5 7.0 -10.0 -16.0];
        
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsToInclude);
        return;
    end

    % Perform the computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic)

        MosaicPoolingOptimizer.performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams, ...
            'maxSFcyclesPerDegree', mosaicParams.subspaceRFmappingParams.maxSFcyclesPerDegree, ...
            'stimSizeDegs', mosaicParams.subspaceRFmappingParams.stimSizeDegs, ....
            'posIncrementDegs', mosaicParams.subspaceRFmappingParams.posIncrementDegs, ...
            'reComputeInputConeMosaicSubspaceRFmappingResponses', ~true, ...
            'reComputeMRGCMosaicSubspaceRFmappingResponses', ~true, ...
            'reComputeRFs', ~true, ...
            'onlyVisualizeOptimallyMappedRFmaps', true);
        
        return;
    end

    % Perform the visualizeVisualRFmapForTargetRGC operation
    if (operationSetToPerformContains.visualizeVisualRFmapForTargetRGC)
        MosaicPoolingOptimizer.performVisualizeVisualRFmapForTargetRGC(mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin);
        return;
    end
end