function run()

    % Ask user which mosaic to use
    [mosaicEcc,  mosaicEccsForSummaryStatistics] = ...
        MosaicPoolingOptimizer.chooseMosaicToUse();

    % Get mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEcc);

    % This controls the visualized spatial support range
    tickSeparationArcMin = 6;
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
        
        xMin = mosaicParams.eccDegs(1) - 1.0;
        xMax = mosaicParams.eccDegs(1) + 1.0;
        yMin = mosaicParams.eccDegs(2) - 0.35;
        yMax = mosaicParams.eccDegs(2) + 0.35;


        identifyPooledCones = true;
        identifyInputCones = ~true;
        plotRFoutlines = true;
        labelRetinalMeridians  = ~true;
        backgroundColor = [1 1 1];
        domainVisualizationLimits = [xMin xMax yMin yMax];
        domainVisualizationTicks = struct('x', -30:0.2:30, 'y', -20:0.2:20);

        if (mosaicParams.eccDegs(1)<0)
            reverseXDir = true;
        else
            reverseXDir = false;
        end

        MosaicPoolingOptimizer.performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones', identifyInputCones, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'plotRFoutlines', plotRFoutlines, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'backgroundColor', backgroundColor, ...
            'reverseXDir', reverseXDir);
        return;
    end

    % Perform the visualizePSFsWithinRGCMosaic operation
    if (operationSetToPerformContains.visualizePSFsWithinRGCMosaic)
        MosaicPoolingOptimizer.performVisualizePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin);
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
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
            'gridlessLineWeightingFuncions', true);

        return;
    end

    % Perfom the summarizeOptimizedSurroundConePoolingModels operation
    if (operationSetToPerformContains.summarizeOptimizedSurroundConePoolingModels)
        MosaicPoolingOptimizer.performSummarizeOptimizedSurroundConePoolingModelsOp(...
            mosaicEccsForSummaryStatistics, ...
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

    % Perform the resetGainOfComputeReadyMRGCMosaic
    if (operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performResetGainOfComputeReadyMRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeConePoolingRFmapAndVisualSTFforTargetRGC operation
    if (operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC)
        
        if (mosaicParams.eccDegs(1)<0)
            reverseXDir = true;
        else
            reverseXDir = false;
        end

        MosaicPoolingOptimizer.performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFuncions', true);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsForSummaryStatistics);
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

   
    if (operationSetToPerformContains.computeVisualRFcenterMapsViaDirectConvolutionWithPSF)
        MosaicPoolingOptimizer.performComputeVisualRFcenterMapsViaDirectConvolutionWithPSF()
    end

end