function [theMosaicFileName, mosaicDirectoryPath] = mosaicFileName(mosaicCenterParams)
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();

    mosaicDirectoryPath = sprintf('%s/productionMidgetRGCMosaics', dropboxDir);
    theMosaicFileName = sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f.mat', ...
        mosaicCenterParams.positionDegs(1), mosaicCenterParams.positionDegs(2), ...
        mosaicCenterParams.sizeDegs(1), mosaicCenterParams.sizeDegs(2));

    theMosaicFileName = fullfile(mosaicDirectoryPath, theMosaicFileName);
end