function mosaicEcc =  chooseMosaicToUse()

    % Select optics to employ
    availableMosaicEccs = {'  0.0', '  2.5', '  7.0', '-10.0', '-20.0'};
    availableMosaicEccsNoSpaces = strrep(availableMosaicEccs, ' ', '');

    theChoice = 'invalid';
    while (~ismember(theChoice , availableMosaicEccsNoSpaces))
       fprintf('\nAvailable mRGC mosaics: \n');
       for iM = 1:numel(availableMosaicEccs)
            fprintf('\t %s\n', availableMosaicEccs{iM});
       end
       theChoice = input('Mosaic to use: ', 's');
    end

    mosaicEcc = str2num(theChoice);
end    