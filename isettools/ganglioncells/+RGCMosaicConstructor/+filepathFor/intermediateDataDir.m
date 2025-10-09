function theIntermediateDataDir = intermediateDataDir()

    p = getpref('isetbio');
    assert(isfield (p, 'rgcResources'), 'Did not find an ''rgcResources'' field in isetbio preferences. Have you run ''RGCMosaicConstructor.helper.utils.generateLocalPrefs() ?''.' );
    assert(isfield (p.rgcResources, 'intermediateDataDir'), 'Did not find an ''intermediateDataDir'' field in isetbio preference ''p.rgcResources''. Have you run ''RGCMosaicConstructor.helper.utils.generateLocalPrefs()?''.' );

    theIntermediateDataDir = p.rgcResources.intermediateDataDir;

    % OLD WAY
    % dropboxDirPath = RGCMosaicConstructor.filepathFor.localDropboxDir();
    % theIntermediateDataDir = fullfile(dropboxDirPath, 'IBIO_rgcMosaicResources', 'ONcenterMidgetRGCmosaics', 'intermediateFiles', 'SLIM');

end