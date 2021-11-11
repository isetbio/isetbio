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

function theOI = wvf2oiSpecial(theWVF, umPerDegree, pupilDiameterMM)

    % Generate oi from the wvf
    theOI = wvf2oi(theWVF);
    
    % Adjust the OI's fNumber and focalLength to be consistent with the
    % micronsPerDegree and pupilDiameter of the WVF
    optics = oiGet(theOI, 'optics');
    focalLengthMM = (umPerDegree * 1e-3) / (2 * tand(0.5));
    focalLengthMeters = focalLengthMM * 1e-3;

    pupilRadiusMeters = (pupilDiameterMM / 2) * 1e-3;
    pupilDiameterMeters = 2 * pupilRadiusMeters;
    optics = opticsSet(optics, 'fnumber', focalLengthMeters / pupilDiameterMeters);
    optics = opticsSet(optics, 'focalLength', focalLengthMeters);
    theOI = oiSet(theOI, 'optics', optics);

    %heightDegs = oiGet(theOI, 'hangular');
    %heightMicrons = oiGet(theOI, 'height')*1e6;
    %fprintf('Achievend microns per deg: %f (desired: %f)\n', heightMicrons/heightDegs, umPerDegree);  
end
