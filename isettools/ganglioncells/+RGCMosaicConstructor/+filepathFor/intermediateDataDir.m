function theIntermediateDataDir = intermediateDataDir()

    p = getpref('isetbio'); 
    theIntermediateDataDir = p.rgcResources.intermediateDataDir;

    % OLD WAY
    % dropboxDirPath = RGCMosaicConstructor.filepathFor.localDropboxDir();
    % theIntermediateDataDir = fullfile(dropboxDirPath, 'IBIO_rgcMosaicResources', 'ONcenterMidgetRGCmosaics', 'intermediateFiles', 'SLIM');

end