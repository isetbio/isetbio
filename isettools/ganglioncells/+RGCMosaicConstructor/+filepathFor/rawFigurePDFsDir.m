% RGCMosaicConstructor.filepathFor.rawFigurePDFsDir
function theRawFiguresDir = rawFigurePDFsDir()

    p = getpref('isetbio'); 
    theRawFiguresDir = p.rgcResources.figurePDFsDir;

    % OLD WAY
    %dropboxDirPath = RGCMosaicConstructor.filepathFor.localDropboxDir();
    %theRawFiguresDir = fullfile(dropboxDirPath, 'ManuscriptSupportMaterials', 'PLOS2024', 'figures', 'SLIM', 'raw');
end