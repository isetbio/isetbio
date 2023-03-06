function inspectVisualSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename)

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');


    % Load the computed mRGC  STF responses 
    load(mRGCMosaicSTFresponsesFilename, 'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    % Find the orientation resulting in the highest extending STF
    % Fit it with a DoG

  

end
