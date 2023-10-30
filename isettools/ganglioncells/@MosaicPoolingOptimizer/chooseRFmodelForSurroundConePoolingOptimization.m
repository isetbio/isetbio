function [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
    chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams, varargin)

    p = inputParser;
    p.addParameter('PackerDaceyH1horizontalCellIndex', 4, @(x)(isscalar(x)||ismember(x,validH1indices)));
    p.parse(varargin{:});
    PackerDaceyH1horizontalCellIndex = p.Results.PackerDaceyH1horizontalCellIndex;

    % Grid sampling 'rectangular' or 'hexagonal'
    gridSamplingScheme = 'rectangular';

    if (mosaicParams.eccDegs(1) == 0)
        gridSamplingScheme = 'hexagonal';
    end

    % RetinalRFmodel params to employ
    retinalRFmodelParams = MosaicPoolingOptimizer.defaultRetinalRFmodelParams;
    retinalRFmodelParams.conePoolingModel = sprintf('arbitraryCenterConeWeights_doubleExpH1cellIndex%dSurroundWeights', PackerDaceyH1horizontalCellIndex);


    [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    optimizedRGCpoolingObjectsFileName = fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName);
end

   