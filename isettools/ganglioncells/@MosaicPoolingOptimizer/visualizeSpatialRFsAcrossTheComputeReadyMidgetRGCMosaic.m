function visualizeSpatialRFsAcrossTheComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilename, pdfFileName)

    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

    targetRGCposition = [5.6 -1.3];
    targetCenterConesNum = 3;
    targetCenterConeMajorityType = cMosaic.LCONE_ID;

    theComputeReadyMRGCmosaic.visualizeSpatialRFnearPosition(...
        targetRGCposition, targetCenterConesNum, ...
        targetCenterConeMajorityType, ...
        'pdfFileName', pdfFileName);

end


    