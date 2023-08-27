function performVisualizePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    [~, ~, pdfsDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Get the optics params
    if (~isempty(theMidgetRGCMosaic.theNativeOpticsParams))
        fprintf('The visualized PSFs correspond for the mosaic''s nativeOpticsParams: \n');
        opticsParams = theMidgetRGCMosaic.theNativeOpticsParams;
    elseif (~isempty(theMidgetRGCMosaic.theCustomOpticsParams))
        fprintf('The visualized PSFs correspond for the mosaic''s customOpticsParams: \n');
        opticsParams = theMidgetRGCMosaic.theCustomOpticsParams;
    else
        fprintf('The visualized PSFs correspond for the mosaic''s defaultOpticsParams: \n');
        opticsParams = theMidgetRGCMosaic.defaultOpticsParams;
    end

    xPosDegs = 0.5*mosaicParams.sizeDegs(1)*[-1 0 1];
    yPosDegs = 0.5*mosaicParams.sizeDegs(2)*[-1 0 1];

    [X,Y] = meshgrid(xPosDegs, yPosDegs);
    eccDegs = [X(:) Y(:)];
    eccDegs = bsxfun(@plus, eccDegs, mosaicParams.eccDegs);

    % Visualize the PSFs
    hFig = theMidgetRGCMosaic.visualizeOpticsAtEccentricities(eccDegs, opticsParams, tickSeparationArcMin);

    % Export fig
    NicePlot.exportFigToPDF(fullfile(pdfsDirectory,'PSFs.pdf'), hFig, 300);

end

