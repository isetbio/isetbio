%RGCMosaicConstructor.helper.utils.feasibilityOfWritingToFile(theInputConeMosaicSTFResponsesFullFileName);

function feasibilityOfWritingToFile(theFileName)

    % Make sure we can write to theFileName
    if (~isfile(theFileName))
        theResult = system(sprintf('touch %s', strrep(theFileName, ' ', '\ ')));
        if (theResult ~= 0)
            error('Cannot write to file %s\n', theFileName);
        end
    end
end


