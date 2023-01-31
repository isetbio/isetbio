function freezeMosaic(mosaicCenterParams, rfModelParams, opticsParams)
    
    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);
    
    % Generate frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);
    
    % Load the center-connected mosaic
    fprintf('\nLoading live mosaic ....\n ');
    load(mosaicFileName, 'theMidgetRGCmosaic');
    fprintf('Done!\n')

    % Freeze it
    fprintf('\nFreezing live mosaic ....\n');
    theMidgetRGCmosaic.freeze();
    fprintf('Done freezing!\n')

    % Save the frozen mosaic
    fprintf('\nSaving frozen mosaic ....\n');
    save(frozenMosaicFileName, 'theMidgetRGCmosaic', '-v7.3');

    % Notify user of the location where the frozen mosaic was saved
    fprintf('Frozen mosaic saved in %s\n', frozenMosaicFileName);
end
