function R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    R2VFTobjFileName = strrep(mosaicFileName, '.mat', sprintf('_R2VFTobjectsForH1cellIndex_%d.mat', H1cellIndex));
end