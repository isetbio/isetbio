function performGenerateCenterConnectedRGCMosaicOp(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', MosaicConnector.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, @isscalar);
    p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
    p.parse(varargin{:});

    % Generate mosaic, its input coneMosaic and connect cones to the RF centers
    theMidgetRGCMosaic = mRGCMosaic(...
        'eccentricityDegs', mosaicParams.eccDegs, ...
        'sizeDegs', mosaicParams.sizeDegs, ...
        'chromaticSpatialVarianceTradeoff', p.Results.chromaticSpatialVarianceTradeoff, ...
        'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', p.Results.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs);

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
