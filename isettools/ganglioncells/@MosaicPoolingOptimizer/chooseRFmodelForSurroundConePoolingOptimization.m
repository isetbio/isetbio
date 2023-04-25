function [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
    chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams)

    % Grid sampling 'rectangular' or 'hexagonal'
    gridSamplingScheme = 'rectangular';

    if (mosaicParams.eccDegs(1) == 0)
        gridSamplingScheme = 'hexagonal';
    end

    % RetinalRFmodel params to employ
    % Change something if we want, like the model name, e.g. choose cell index 3,
    % 'arbitraryCenterConeWeights_doubleExpH1cellIndex3SurroundWeights', ... 
    H1cellIndex = 4;
    retinalRFmodelParams = MosaicPoolingOptimizer.defaultRetinalRFmodelParams;
    retinalRFmodelParams.conePoolingModel = sprintf('arbitraryCenterConeWeights_doubleExpH1cellIndex%dSurroundWeights', H1cellIndex);


    [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    optimizedRGCpoolingObjectsFileName = fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName);
end

   