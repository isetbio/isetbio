function [theMosaicFileName, mosaicDirectoryPath] = mosaicFileName(mosaicParams)
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
    mosaicDirectoryPath = sprintf('%s/productionMidgetRGCMosaics', dropboxDir);
    theMosaicFileName = sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f.mat', ...
        mosaicParams.positionDegs(1), mosaicParams.positionDegs(2), ...
        mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2));

    theMosaicFileName = fullfile(mosaicDirectoryPath, theMosaicFileName);
end