%
% RGCMosaicConstructor.filepathFor.rawFigurePDFsDir
%
function theRawFiguresDir = rawFigurePDFsDir()

    rgcResources = RGCMosaicConstructor.helper.utils.rgcResources();
    errorMessage = sprintf('\nDid not find a ''figurePDFsDir'' field in isetbio preference ''p.rgcResources''.\nHave you run ''RGCMosaicConstructor.helper.utils.generateLocalPrefs()?''.');
    assert(isfield (rgcResources, 'figurePDFsDir'), errorMessage);
    theRawFiguresDir = rgcResources.figurePDFsDir;

 end