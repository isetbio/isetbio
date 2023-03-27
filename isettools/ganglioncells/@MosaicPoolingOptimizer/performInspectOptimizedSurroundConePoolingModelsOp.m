function performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');


    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index was used to optimize the RF
    % surround pooling model models
    [~, ~, optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);


    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'samplingScheme', 'rectangular', ...
         'generateSamplingGrids', true, ...
         'visualizeSamplingGrids', false);

    gridNodesToInspect = input('Enter grid node to inspect. Hit enter to inspect all.: ');
    if (isempty(gridNodesToInspect))
       gridNodesToInspect = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum;
    end

    tickSeparationArcMin = 6;
    visualizedSpatialFrequencyRange = []; %[0.1 30];

    for iNode = 1:numel(gridNodesToInspect)
       gridNodeIndex = gridNodesToInspect(iNode);
       theMosaicPoolingOptimizerOBJ.inspect(...
           gridNodeIndex, ...
           opticsParams, ...
           fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
           'tickSeparationArcMin', tickSeparationArcMin, ...
           'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange ...
           );
    end

end
