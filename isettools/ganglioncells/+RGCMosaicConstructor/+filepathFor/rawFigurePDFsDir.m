function theRawFiguresDir = rawFigurePDFsDir()
    dropboxDirPath = RGCMosaicConstructor.filepathFor.localDropboxDir();
    theRawFiguresDir = fullfile(dropboxDirPath, 'ManuscriptSupportMaterials', 'PLOS2024', 'figures', 'SLIM', 'raw');
end