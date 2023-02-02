function generateCenterSurroundConePoolingKernels(mosaicCenterParams,  rfModelParams, opticsParams)

    midgetRGCMosaicInspector.say('Generating center-surround cone pooling kernels for all RGCs in the mosaic');

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);
    
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(...
        mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the computed R2VFTobjects
    load(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

    % Generate the center/surround RFs
    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRFs(theRTFVTobjList, ...
                        theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
                        theVisualSTFSurroundToCenterRcRatioGrid, theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid);

    % Save the updated midgetRGCmosaic which now includes  the computed
    % RTVFTobjList as well as the different grids:
    %  - 'theOpticsPositionGrid'
    %  - 'theConesNumPooledByTheRFcenterGrid'
    %  - 'theVisualSTFSurroundToCenterRcRatioGrid'
    %  - 'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid'
    save(mosaicFileName, 'theMidgetRGCmosaic', '-v7.3');

    fprintf('MidgetRGCMosaic with center/surround pooling weights and RTVFTobjList saved in: %s\n', mosaicFileName);

end
