function theIntermediateDataDir = intermediateDataDir()
    dropboxDirPath = RGCMosaicConstructor.filepathFor.localDropboxDir();
    theIntermediateDataDir = fullfile(dropboxDirPath, 'IBIO_rgcMosaicResources', 'ONcenterMidgetRGCmosaics', 'intermediateFiles');
end