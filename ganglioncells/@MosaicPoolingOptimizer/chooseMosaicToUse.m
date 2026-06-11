function [mosaicHorizontalEccentricityDegs, mosaicEccsForSummaryStatistics] =  chooseMosaicToUse()

    % Select optics to employ
    availableMosaicEccs = {'  0.0', '  2.5', '  7.0', '-10.0', '-16.0'};
    for iEcc = 1:numel(availableMosaicEccs)
        mosaicEccsForSummaryStatistics(iEcc) = str2num(availableMosaicEccs{iEcc});
    end

    availableMosaicEccsNoSpaces = strrep(availableMosaicEccs, ' ', '');

    theChoice = 'invalid';
    while (~ismember(theChoice , availableMosaicEccsNoSpaces))
       fprintf('\nAvailable mRGC mosaics: \n');
       for iM = 1:numel(availableMosaicEccs)
            fprintf('\t %s\n', availableMosaicEccs{iM});
       end
       theChoice = input('Mosaic to use: ', 's');
    end

    mosaicHorizontalEccentricityDegs = str2num(theChoice);
end    