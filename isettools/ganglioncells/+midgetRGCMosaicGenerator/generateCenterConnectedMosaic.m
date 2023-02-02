function generateCenterConnectedMosaic(mosaicCenterParams)

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);
    
    fprintf('Mosaic will be saved in: %s', mosaicFileName);
 
    midgetRGCMosaicInspector.say('Generating center-connected midgetRGCMosaic. Please wait');

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