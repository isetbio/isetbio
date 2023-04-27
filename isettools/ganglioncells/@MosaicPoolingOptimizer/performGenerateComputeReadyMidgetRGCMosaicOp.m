function performGenerateComputeReadyMidgetRGCMosaicOp(mosaicParams)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate and set the optics
    theMidgetRGCMosaic.setTheOptics(opticsParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);


    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);
    
    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'samplingScheme', gridSamplingScheme, ...
         'generateSamplingGrids', true, ...
         'visualizeSamplingGrids', false);


    % Ask user whether to visualize the interpolation process
    visualizeInterpolation = input('Visualize interpolation? [y = YES] : ', 's');
    if (strcmp(visualizeInterpolation, 'y'))
        visualizeInterpolation = true;
    else
        visualizeInterpolation = false;
    end
        
    theMosaicPoolingOptimizerOBJ.generateComputeReadyMidgetRGCMosaic(...
            optimizedRGCpoolingObjectsFileName, ...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            visualizeInterpolation);

end
