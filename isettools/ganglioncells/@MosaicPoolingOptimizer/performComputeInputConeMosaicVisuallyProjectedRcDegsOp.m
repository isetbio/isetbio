function performComputeInputConeMosaicVisuallyProjectedRcDegsOp(mosaicParams)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Ask the user about what optics to use for computing the input cone
    % mosaic STF responses
    [~, opticsToEmploy, coneMosaicSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic and no sampling grid
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'generateSamplingGrids', false);


    % Compute the input cone mosaic visually projectedRcDegs by fitting a Gaussian to the
    % previously computed input cone mosaic visual STFs
    theMosaicPoolingOptimizerOBJ.computeInputConeMosaicVisuallyProjectedRcDegs(...
            coneMosaicSTFresponsesFileName, ...
            'useParfor', true, ...
            'visualizedResponses', ~true, ...
            'opticsToEmploy', opticsToEmploy);

end
