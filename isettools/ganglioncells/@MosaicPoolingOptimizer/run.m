function run()
    
    % Eccentricities of optimized mosaics
    mosaicEcc = 2.5;
    mosaicEcc = 7.0;
    mosaicEcc = -10.0;

    % Get mosaic params
    mosaicParams = getMosaicParams(mosaicEcc);

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
        MosaicPoolingOptimizer.performOptimizeSurroundConePoolingModelsOp(mosaicParams);
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
    if (operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic )
        MosaicPoolingOptimizer.performComputeVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

end


function mosaicParams = getMosaicParams(mosaicEcc)
    switch (mosaicEcc)
        case 0
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
            % which covers the [1 - 4] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [0 0], ...
                'sizeDegs', [3 3]);

        case 2.5
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
            % which covers the [1 - 4] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [2.5 0], ...
                'sizeDegs', [3 3]);

        case 7.0
            % Mosaic params to employ. This is for the 7.0 deg - centered mosaic
            % which covers the [4-10] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [7 0], ...
                'sizeDegs', [6 3]);

        case -10.0
            mosaicParams = struct(...
                'eccDegs', [-10 0], ...
                'sizeDegs', [6 3]);

        otherwise
            error('No data for this eccentricity')
    end
end
