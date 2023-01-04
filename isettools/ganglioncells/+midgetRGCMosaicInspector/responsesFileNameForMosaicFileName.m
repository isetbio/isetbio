function  responsesFileName = responsesFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    responsesFileName = strrep(mosaicFileName, '.mat', sprintf('_H1cellIndex%d_Responses.mat', H1cellIndex));
end
