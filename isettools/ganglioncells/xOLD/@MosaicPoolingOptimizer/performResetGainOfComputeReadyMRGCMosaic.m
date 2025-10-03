function performResetGainOfComputeReadyMRGCMosaic(mosaicParams)

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    computeReadyMosaicFileName = fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName);
    load(computeReadyMosaicFileName, 'theComputeReadyMRGCmosaic');

    % Set the peak gains based on the net cone weights to the cennter(default)
    method = '1/integrated center cone weights';
    gain = 1.0;
    methodParams = gain;
    theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);

    % Save the updated mosaic
    save(computeReadyMosaicFileName, 'theComputeReadyMRGCmosaic', '-v7.3');
    fprintf('The default rgcRFgains - compute-ready mRGCMosaic was exported in %s.\n', computeReadyMosaicFileName);
end