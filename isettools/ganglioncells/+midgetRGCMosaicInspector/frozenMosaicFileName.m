function theFrozenMosaicFileName = frozenMosaicFileName(mosaicCenterParams, H1cellIndex)
    theLiveMosaicFileName = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);
    theFrozenMosaicFileName = strrep(theLiveMosaicFileName, '.mat', sprintf('_H1cellIndex_%d_Frozen.mat', H1cellIndex));
    % Separate directory
    theFrozenMosaicFileName = strrep(theFrozenMosaicFileName, 'MRGCmosaic', 'frozenMosaics/MRGCmosaic');
end