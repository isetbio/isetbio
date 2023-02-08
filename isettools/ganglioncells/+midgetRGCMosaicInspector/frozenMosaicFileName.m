function theFrozenMosaicFileName = frozenMosaicFileName(mosaicCenterParams, H1cellIndex, opticsParams)

    opticsPostFix = sprintf('%s_Rank_%d_Pupil_%1.1f', ...
                    opticsParams.ZernikeDataBase, opticsParams.subjectRankOrder, opticsParams.pupilDiameterMM);

    modelEncodingPostFix = sprintf('_%s_H1cellIndex_%d', opticsPostFix, H1cellIndex);

    theLiveMosaicFileName = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);
    theFrozenMosaicFileName = strrep(theLiveMosaicFileName, '.mat', sprintf('%s_Frozen.mat', modelEncodingPostFix));
    
    % Separate directory
    theFrozenMosaicFileName = strrep(theFrozenMosaicFileName, 'MRGCmosaic', 'frozenMosaics/MRGCmosaic');
end