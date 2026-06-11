function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = ...
    diffractionLimitedOptics(pupilDiameterMM, wavelengthsListToCompute, imagingPeakWavelength, ...
    micronsPerDegree, opticalDefocusDiopters, varargin)

    p = inputParser;
    p.addParameter('noLCA', false, @islogical);
    p.parse(varargin{:});
    noLCA = p.Results.noLCA;
    
    wavefrontSpatialSamples = 401;
    zCoeffs = zeros(1,22);

    % Add custom optical defocus
    if (opticalDefocusDiopters ~= 0)
        refractiveErrorMicrons = wvfDefocusDioptersToMicrons(opticalDefocusDiopters, pupilDiameterMM);
        zCoeffs(1+4) = refractiveErrorMicrons;
    end
    
    % Compute the optics
    [thePSF, ~, ~,~, psfSupportMinutesX, psfSupportMinutesY, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             pupilDiameterMM, ...
             pupilDiameterMM, imagingPeakWavelength, false, ...
             'doNotZeroCenterPSF', ~false, ...
             'micronsPerDegree', micronsPerDegree, ...
             'upsampleFactor', 2, ...
             'name', 'DiffractionLimitedOptics');
         
    % Remove wavelength-dependent defocus if noLCA is set
    if (noLCA)
        % Set all PSFs to the PSF at the in-focus wavelenth
        [~,wTarget] = min(abs(wavelengthsListToCompute-imagingPeakWavelength));
        targetPSF = thePSF(:,:,wTarget);
        for waveIndex = 1:numel(wavelengthsListToCompute)
            theWVF.psf{waveIndex} = targetPSF;
            thePSF(:,:,waveIndex) = targetPSF;
        end
    end
    
    % Generate the OI from the wavefront map
    theOI = wvf2oiSpecial(theWVF, micronsPerDegree, pupilDiameterMM);
    
    psfSupportWavelength = wavelengthsListToCompute;
end
