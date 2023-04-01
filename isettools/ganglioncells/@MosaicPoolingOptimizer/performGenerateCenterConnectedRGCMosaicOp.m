function performGenerateCenterConnectedRGCMosaicOp(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', MosaicConnector.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, @isscalar);
    p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
    p.parse(varargin{:});

    % Retrieve input cone mosaic
    existingMRGCMosaicFileName = fullfile(MosaicPoolingOptimizer.localDropboxPath(), 'productionMidgetRGCMosaics/MosaicOptimizerResources/mRGCMosaicEcDegs(-20.0_0.0)_SizeDegs(10.0_5.0).mat');
    load(existingMRGCMosaicFileName, 'theMidgetRGCMosaic');
    theInputConeMosaic = theMidgetRGCMosaic.inputConeMosaic;

    % Generate mosaic, its input coneMosaic and connect cones to the RF centers
    theMidgetRGCMosaic = mRGCMosaic(...
        'eccentricityDegs', mosaicParams.eccDegs, ...
        'sizeDegs', mosaicParams.sizeDegs, ...
        'withInputConeMosaic', theInputConeMosaic, ...
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
