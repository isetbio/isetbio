function [theLconeWeightedPSFstruct, theMconeWeightedPSFstruct, theSconeWeightedPSFstruct] = coneWeightedPSFs(thePSF)
    % Load the 2-deg Stockman cone fundamentals on wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), thePSF.supportWavelength);
    theLconeWeightingSpectrum = squeeze(coneFundamentals(:,1));
    theMconeWeightingSpectrum = squeeze(coneFundamentals(:,2));
    theSconeWeightingSpectrum = squeeze(coneFundamentals(:,3));
    theLconeWeightingSpectrum = theLconeWeightingSpectrum / sum(theLconeWeightingSpectrum(:));
    theMconeWeightingSpectrum = theMconeWeightingSpectrum / sum(theMconeWeightingSpectrum(:));
    theSconeWeightingSpectrum = theSconeWeightingSpectrum / sum(theSconeWeightingSpectrum(:));

    nRows = size(thePSF.data,1);
    nCols = size(thePSF.data,2);
    theLconeWeightedPSF = zeros(nRows, nCols);
    theMconeWeightedPSF = zeros(nRows, nCols);
    theSconeWeightedPSF = zeros(nRows, nCols);
    for iWave = 1:numel(thePSF.supportWavelength)
        thePSFslice = squeeze(thePSF.data(:,:,iWave));
        theLconeWeightedPSF = theLconeWeightedPSF + theLconeWeightingSpectrum(iWave) * thePSFslice;
        theMconeWeightedPSF = theMconeWeightedPSF + theMconeWeightingSpectrum(iWave) * thePSFslice;
        theSconeWeightedPSF = theSconeWeightedPSF + theSconeWeightingSpectrum(iWave) * thePSFslice;
    end
    maxPSF = max([max(theLconeWeightedPSF(:)) max(theMconeWeightedPSF(:)) max(theSconeWeightedPSF(:))]);

    theLconeWeightedPSFstruct.data = reshape(theLconeWeightedPSF, [nRows nCols 1]) / maxPSF;
    theMconeWeightedPSFstruct.data = reshape(theMconeWeightedPSF, [nRows nCols 1]) / maxPSF;
    theSconeWeightedPSFstruct.data = reshape(theSconeWeightedPSF, [nRows nCols 1]) / maxPSF;
    theLconeWeightedPSFstruct.supportWavelength = 1;
    theMconeWeightedPSFstruct.supportWavelength = 1;
    theSconeWeightedPSFstruct.supportWavelength = 1;

    theLconeWeightedPSFstruct.supportX = thePSF.supportX;
    theLconeWeightedPSFstruct.supportY = thePSF.supportY;
    theMconeWeightedPSFstruct.supportX = thePSF.supportX;
    theMconeWeightedPSFstruct.supportY = thePSF.supportY;
    theSconeWeightedPSFstruct.supportX = thePSF.supportX;
    theSconeWeightedPSFstruct.supportY = thePSF.supportY;
end