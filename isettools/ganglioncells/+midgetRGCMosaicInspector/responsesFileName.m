function  responsesFileName = responsesFileName(mosaicFileName, opticsPositionDegs)
    responsesFileNamePostFix = sprintf('_ResponsesForOpticsAt_%2.2f_%2.2fdegs.mat', ...
        opticsPositionDegs(1), opticsPositionDegs(2));
    
    responsesFileName = strrep(mosaicFileName, '.mat', responsesFileNamePostFix);
end
