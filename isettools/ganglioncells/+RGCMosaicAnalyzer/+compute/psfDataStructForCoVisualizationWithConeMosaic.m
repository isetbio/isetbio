%
% RGCMosaicAnalyzer.compute.psfDataStructForCoVisualizationWithConeMosaic
%

function psfAndMtfDataStruct = psfDataStructForCoVisualizationWithConeMosaic(theOI, domain, theVisualizedWavelength)

    % Get the optics data
    optics = oiGet(theOI, 'optics');

    % Get the focal length
    focalLengthMeters = opticsGet(optics, 'focal length');
    micronsPerDegree = focalLengthMeters*tand(1)*1e6;

    % Get the wavelength support
    wavelengthSupport = opticsGet(optics, 'otf wave');
    [~,theVisualizedWavelengthIndex] = min(abs(wavelengthSupport-theVisualizedWavelength));

    % Get the psf data
    psf = opticsGet(optics, 'psf data');
    theTargetWavelengthPSF = squeeze(psf.psf(:,:,theVisualizedWavelengthIndex));

    % Compute the vLambda weighted PSF
    psfDataStruct.data = psf.psf;
    psfDataStruct.supportWavelength = wavelengthSupport;
    theVlamdaWeightedPSF = RGCMosaicAnalyzer.compute.vLambdaWeightedPSF(psfDataStruct);

    % Get the psf support
    psfSupportXmicrons = squeeze(psf.xy(1,:,1));
    psfSupportYmicrons = squeeze(psf.xy(:,2,2));
   
    % Get the OTF
    theOTF = opticsGet(optics,'otf');
    theVisualizedOTF = squeeze(theOTF(:,:,theVisualizedWavelengthIndex));
    theVisualizedMTF = fftshift(abs(theVisualizedOTF));

    % Circularly symmetric MTF
    theVisualizedCircularlySymmetricMTF = psfCircularlyAverage(theVisualizedMTF);

    % Get the OTF support
    otfSupportCyclesPerMM = opticsGet(optics,'otf fx', 'mm');
    otfSupportCyclesPerMicron = otfSupportCyclesPerMM * 1e-3;
    otfSupportCyclesPerDeg = otfSupportCyclesPerMicron * micronsPerDegree;

    if (strcmp(domain, 'microns'))
        psfAndMtfDataStruct.supportXmicrons = psfSupportXmicrons;
        psfAndMtfDataStruct.supportYmicrons = psfSupportYmicrons;
    else
        psfAndMtfDataStruct.supportXdegs = psfSupportXmicrons / micronsPerDegree;
        psfAndMtfDataStruct.supportYdegs = psfSupportYmicrons / micronsPerDegree;
    end

    psfAndMtfDataStruct.psfSliceAtTargetWavelength = theTargetWavelengthPSF;
    psfAndMtfDataStruct.psfVlambdaWeighted = theVlamdaWeightedPSF;
    psfAndMtfDataStruct.mtf = theVisualizedMTF;
    psfAndMtfDataStruct.mtfCircularlySymmetric = theVisualizedCircularlySymmetricMTF;
    psfAndMtfDataStruct.mtfSupportCyclesPerDeg = otfSupportCyclesPerDeg;
end


