function theFileName = encodeMSequenceSpecificInfoToFilename(theFileName, ...
        ternaryInsteadOfBinaryMsequence, mSequenceBitLength, rfPixelsAcross)

    % Encode m-sequence type in responses filename
    if (ternaryInsteadOfBinaryMsequence)
        theFileName = strrep(theFileName, 'MSequence', sprintf('Ternary%2.0fBitMSequence', mSequenceBitLength));
    else
        theFileName = strrep(theFileName, 'MSequence', sprintf('Binary%2.0fBitMSequence', mSequenceBitLength));
    end

    % Encode rfPixelsAcross
    theFileName = strrep(theFileName, 'MSequence', sprintf('MSequence%dx%dPixels', rfPixelsAcross, rfPixelsAcross));
end
