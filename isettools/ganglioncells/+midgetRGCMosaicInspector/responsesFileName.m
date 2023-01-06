function  responsesFileName = responsesFileName(mosaicFileName, H1cellIndex, opticsPositionDegs)
    responsesFileName = strrep(mosaicFileName, '.mat', sprintf('_H1cellIndex%d_ResponsesForOpticsAt_2.2f_%2.2degs.mat', ...
        H1cellIndex, opticsPositionDegs(1), opticsPositionDegs(2)));
end
