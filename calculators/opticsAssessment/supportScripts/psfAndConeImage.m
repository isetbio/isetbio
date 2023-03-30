function [psfImage, coneImage, mosaicNyquistFrequencyCPD, psfSupportArcMin, zCoeffs] = psfAndConeImage(mosaicEcc, opticsParams)

    mosaicEccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(mosaicEcc);
    coneSpacingDegs = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(...
        opticsParams.whichEye, mosaicEccMicrons, 'cones', false);
    
    mosaicSizeDegs = coneSpacingDegs * 10;
    wavelength = 550;
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', mosaicSizeDegs*[1 1], ...    % SIZE in degs
        'eccentricityDegs', mosaicEcc, ...       % ECC in degs
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', wavelength, ...
        'opticalImagePositionDegs', 'mosaic-centered', ...
        'whichEye', opticsParams.whichEye ...
        );
    
    if (isempty(cm.coneRFspacingsDegs))
        coneImage = [];
        psfImage = [];
        mosaicNyquistFrequencyCPD = [];
        psfSupportArcMin = []; 
        zCoeffs = [];
        return;
    end
    
    
    % Compte Nyquist frequency based on mean cone spacing in the mosaic
    mosaicNyquistFrequencyCPD = 1/(2*mean(cm.coneRFspacingsDegs));

    wavefrontSpatialSamples = 901;
    
    switch (opticsParams.zernikeDataBase)
        case 'Artal2012'
        [oiEnsemble, psfEnsemble, zCoeffs] = ...
                cm.oiEnsembleGenerate(mosaicEcc, ...
                'zernikeDataBase', opticsParams.zernikeDataBase, ...
                'subjectID', opticsParams.subjectID, ...
                'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                'flipPSFUpsideDown', opticsParams.flipPSFUpsideDown);
        psf = psfEnsemble{1};
                

        case 'Polans2015'
        [oiEnsemble, psfEnsemble, zCoeffs] = ...
                cm.oiEnsembleGenerate(mosaicEcc, ...
                'zernikeDataBase', opticsParams.zernikeDataBase, ...
                'subjectID', opticsParams.subjectID, ...
                'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                'flipPSFUpsideDown', opticsParams.flipPSFUpsideDown);
        psf = psfEnsemble{1};
        
        case 'Thibos2002'
        [oiEnsemble, psfEnsemble, zCoeffs] = ...
                cm.oiEnsembleGenerate(mosaicEcc, ...
                'zernikeDataBase', opticsParams.zernikeDataBase, ...
                'subjectID', opticsParams.subjectID, ...
                'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                'flipPSFUpsideDown', opticsParams.flipPSFUpsideDown);
        psf = psfEnsemble{1};
        
        otherwise
            error('Unkniwn database: ''%s''.', opticsParams.zernikeDataBase);
    end
    
    if isnan(min(psf.data(:))) || isnan(max(psf.data(:)))
        psfImage = [];
        coneImage = []; 
        psfSupportArcMin = [];
        return;
    end
    
    % Make PSF image
    [~,idx] = min(abs(psf.supportWavelength-wavelength));
    psfImage = squeeze(psf.data(:,:,idx));
     
    psfSupportMicrons = cm.micronsPerDegree * psf.supportX/60;
    [Xorig,Yorig] = meshgrid(psfSupportMicrons,psfSupportMicrons);
     
    interpolationSamples = 1024;
    psfSupportMicrons = linspace(psfSupportMicrons(1), psfSupportMicrons(end), interpolationSamples);
    [Xhires,Yhires] = meshgrid(psfSupportMicrons,psfSupportMicrons);
    psfImage = interp2(Xorig, Yorig, psfImage, Xhires, Yhires);
    psfImage = psfImage / max(psfImage(:));
    psfSupportArcMin = psfSupportMicrons / cm.micronsPerDegree * 60;
    
    % Make cone image
    coneImage = Xhires*0;
    r = sqrt(Xhires.^2 + Yhires.^2);
    
    % Cone aperture radius
    apertureRadiusMicrons = 0.5*cm.coneApertureToDiameterRatio*mean(cm.coneRFspacingsMicrons);
    coneImage(find(r<=apertureRadiusMicrons)) = 1;
end