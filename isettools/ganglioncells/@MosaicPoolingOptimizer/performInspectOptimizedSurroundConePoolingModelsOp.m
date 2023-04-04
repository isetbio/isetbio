function performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams, varargin)

   % Parse optional input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
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
    [~, gridSamplingScheme,  optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);


    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'samplingScheme', gridSamplingScheme, ...
         'generateSamplingGrids', true, ...
         'visualizeSamplingGrids', false);

    queryString = sprintf('\nEnter grid node to inspect [%d-%d]. Alternatively, hit enter to inspect multiple nodes: ', ...
        1, theMosaicPoolingOptimizerOBJ.gridNodesNum);
    gridNodesToInspect = input(queryString, 's');
    if (isempty(gridNodesToInspect))
        gridNodesToInspect = MosaicPoolingOptimizer.gridNodesToOptimize();
        if (isempty(gridNodesToInspect))
             gridNodesToInspect = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum;
        end

    else
       gridNodesToInspect = str2num(gridNodesToInspect);
    end


    for iNode = 1:numel(gridNodesToInspect)

       if (iscell(gridNodesToInspect))
            gridNodeIndex = gridNodesToInspect{iNode}.number;
       else
            gridNodeIndex = gridNodesToInspect(iNode);
       end

       theMosaicPoolingOptimizerOBJ.inspect(...
           gridNodeIndex, ...
           opticsParams, ...
           fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
           'tickSeparationArcMin', tickSeparationArcMin, ...
           'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange ...
           );
    end

end
