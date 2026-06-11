function [theComputeReadyMRGCmosaic, fullPathToComputeReadyMRGCMosaic] = ...
    queryUserAndLoadComputeReadyMRGCMosaic(mosaicParams)

     % Ask the user what optics were used for computing the compute-ready MRGC mosaic
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParamsForComputeReadyMosaic = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index was used for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(...
        mosaicParams, opticsParamsForComputeReadyMosaic);

    % Generate the filename of the compute-ready mRGCMosaic to use
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParamsForComputeReadyMosaic, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready mRGCMosaic
    fullPathToComputeReadyMRGCMosaic = fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName);
    load(fullPathToComputeReadyMRGCMosaic, 'theComputeReadyMRGCmosaic');

end