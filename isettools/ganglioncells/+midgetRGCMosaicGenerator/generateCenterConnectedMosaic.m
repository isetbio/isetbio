function generateCenterConnectedMosaic(mosaicCenterParams)

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.generateMosaicFileName(mosaicCenterParams);
    
    fprintf('Mosaic will be saved in: %s', mosaicFileName);
    fprintf('Generating mosaic. Please wait ...');
    
    % Generate mRGC mosaic
    theMidgetRGCmosaic = midgetRGCMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'whichEye', mosaicCenterParams.whichEye, ...
                'eccentricityDegs', mosaicCenterParams.positionDegs, ...
                'sizeDegs', mosaicCenterParams.sizeDegs ...
                );
    % Save the center-connected mosaic
    save(mosaicFileName, 'theMidgetRGCmosaic', '-v7.3');
end