%
% RGCMosaicConstructor.filepathFor.rawFigurePDFsDir
%
function theRawFiguresDir = rawFigurePDFsDir()

    rgcResources = RGCMosaicConstructor.helper.utils.rgcResources();
    assert(isfield (rgcResources, 'figurePDFsDir'), 'Did not find a ''figurePDFsDir'' field in isetbio preference ''p.rgcResources''. Have you run ''RGCMosaicConstructor.helper.utils.generateLocalPrefs()?''.' );
    theRawFiguresDir = rgcResources.figurePDFsDir;

 end