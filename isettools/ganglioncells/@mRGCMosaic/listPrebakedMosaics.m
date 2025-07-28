function [themRGCmosaicFileNames, prebakedMRGCMosaicDir] = listPrebakedMosaics()
    
    prebakedMRGCMosaicDir = 'isettools/ganglioncells/data/prebakedRGCmosaics/ONmRGCmosaics';
    prebakedMRGCMosaicDir = fullfile(isetbioRootPath, prebakedMRGCMosaicDir);
    listing = dir(prebakedMRGCMosaicDir);

    folders = struct2table(listing);
    themRGCmosaicFileNames = folders(~matches(folders.name,[".",".."]),:)
end