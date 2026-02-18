%
% RGCMosaicConstructor.filepathFor.intermediateDataDir
% 
function theIntermediateDataDir = intermediateDataDir()

    rgcResources = RGCMosaicConstructor.helper.utils.rgcResources();
    assert(isfield(rgcResources, 'intermediateDataDir'), 'Did not find a ''intermediateDataDir'' field in isetbio preference ''p.rgcResources''. Have you run ''RGCMosaicConstructor.helper.utils.generateLocalPrefs()?''.' );
    theIntermediateDataDir = rgcResources.intermediateDataDir;

end