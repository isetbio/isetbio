function theMRGCMosaic = loadPreComputedMRGCMosaic(horizontalEccDegs)

    % Choose one of the computed mRGCMosaics
    % 1. mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(horizontalEccDegs);

    % 2. optics params
    opticsParams = MosaicPoolingOptimizer.getOpticsParams(mosaicParams);

    % 3. retinal RF pooling params
    retinalRFmodelParams = MosaicPoolingOptimizer.getSurroundParams(mosaicParams, opticsParams);

    % Load the mosaic
    theMRGCMosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
        mosaicParams, opticsParams, retinalRFmodelParams);
end