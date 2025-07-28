function theComputeReadyMRGCmosaic = loadComputeReadyRGCMosaic(mosaicParams, opticsParams, retinalRFmodelParams)

    % Load the compute-ready RGC mosaic
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);
    
    sourceRGCMosaicFileName = fullfile(resourcesDirectory, computeReadyMosaicFileName);

    fprintf('Will try to load %s ... \n', sourceRGCMosaicFileName)
    fprintf('from %s ... \n', resourcesDirectory);

    if ((~isfolder(resourcesDirectory)) || (~isfile(sourceRGCMosaicFileName)))
        % Mosaic was not found. List available compute-ready mosaics of the desired type
        mRGCMosaic.availableComputeReadyMosaics(mosaicParams.rgcType);
    end

    % Check that the mosaic directory exists
    assert(isfolder(resourcesDirectory), sprintf('Mosaic directory (''%s'') not found.', resourcesDirectory));

    % Check that the mosaic file exists
    assert(isfile(sourceRGCMosaicFileName), sprintf('Mosaic file (''%s'') not found.', sourceRGCMosaicFileName));

    % Mosaic file found, so load the data
    load(sourceRGCMosaicFileName, 'theComputeReadyMRGCmosaic');
    fprintf('Compute-ready mRGCmosaic loaded!\n');
end