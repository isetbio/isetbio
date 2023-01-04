function generateCenterSurroundConePoolingKernels(mosaicCenterParams, mosaicSurroundParams)

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.generateMosaicFileName(...
        mosaicCenterParams);
    

    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileNameForMosaicFileName(...
        mosaicFileName, mosaicSurroundParams.H1cellIndex);

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

end
