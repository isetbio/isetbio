function frozenMosaicFileName = frozenMosaicFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    frozenMosaicFileName = strrep(mosaicFileName, '.mat', sprintf('_H1cellIndex_%d_Frozen.mat', H1cellIndex));
    % Separate directory
    frozenMosaicFileName = strrep(frozenMosaicFileName, 'MRGCmosaic', 'frozenMosaics/MRGCmosaic');
end