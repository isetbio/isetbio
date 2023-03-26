function generateCenterConnectedRGCMosaicOperation(mosaicParams)
    % Generate mosaic, its input coneMosaic and connect cones to the RF centers
    theMidgetRGCMosaic = mRGCMosaic(...
        'eccentricityDegs', mosaicParams.eccDegs, ...
        'sizeDegs', mosaicParams.sizeDegs);

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Save the generated center-only connected mRGCmosaic
    save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Visualize input cone mosaic
    theMidgetRGCMosaic.inputConeMosaic.visualize()

    % Visualize cone pooling by the RF centers
    theMidgetRGCMosaic.visualize();
end
