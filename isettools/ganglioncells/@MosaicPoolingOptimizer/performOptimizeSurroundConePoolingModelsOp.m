function performOptimizeSurroundConePoolingModelsOp(mosaicParams)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    
    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    [opticsParams, ~, coneMosaicSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Select the nodes at which to optimize the surround pooling
    gridNodesToOptimize = MosaicPoolingOptimizer.gridNodesToOptimize();

    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'samplingScheme', gridSamplingScheme, ...
         'generateSamplingGrids', true, ...
         'visualizeSamplingGrids', true);

    if (~isempty(gridNodesToOptimize))&&(iscell(gridNodesToOptimize))
       % Optimizing an arbitrary selection of grid nodes
       allGridNodesToCompute = gridNodesToOptimize;
    else
       % Optimizing all grid nodes 
       allGridNodesToCompute = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum;
    end

    % Query user which node indices to compute
    nodeIndicesToWorkOn = input('Which node indices set to work on ? (even, odd, or all) : ', 's');
    switch nodeIndicesToWorkOn
        case 'even'
            gridNodesToCompute = allGridNodesToCompute(2:2:numel(allGridNodesToCompute));
        case 'odd'
            gridNodesToCompute = allGridNodesToCompute(1:2:numel(allGridNodesToCompute));
        case 'all'
            gridNodesToCompute = allGridNodesToCompute(1:1:numel(allGridNodesToCompute));
        otherwise
            error('grid nodes must be either ''even'' or ''odd''.');
    end

    % Optimize surround cone pooling weights for each grid node 
    for iNode = 1:numel(gridNodesToCompute)

        if (iscell(gridNodesToCompute))
            gridNodeIndex = gridNodesToCompute{iNode}.number;
            whichConeType = gridNodesToCompute{iNode}.coneType;
        else
            gridNodeIndex = gridNodesToCompute(iNode);
            whichConeType = [cMosaic.LCONE_ID cMosaic.MCONE_ID];
        end
        
        % Fitting options
        multiStartsNumRetinalPooling = 12;
        multiStartsNumDoGFit = 128;

        % More weight for matching the Rs/Rc ratio
        rmseWeightForRsRcResidual = 2.0;
        rmseWeightForSCintSensResidual = 1.0;

        theMosaicPoolingOptimizerOBJ.compute(gridNodeIndex, whichConeType, opticsParams, ...
            fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
            fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
            'multiStartsNumDoGFit', multiStartsNumDoGFit, ...
            'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
            'rmseWeightForRsRcResidual', rmseWeightForRsRcResidual, ...
            'rmseWeightForSCintSensResidual', rmseWeightForSCintSensResidual, ...
            'retinalRFmodelParams', retinalRFmodelParams, ...
            'displayFittingProgress', true);
    end

end
