function performVisualizeVisualRFmapForTargetRGC(mosaicParams)

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

     % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    load(fullfile(resourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

    % Generate filename for the mRGCMosaic subspace responses
    [mRGCMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    % PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory');

    MosaicPoolingOptimizer.visualizeVisualRFmapsForMultipleTargetRGCs(...
            theComputeReadyMRGCmosaic, ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(pdfDirectory, 'visualRFmap.pdf'));
end
