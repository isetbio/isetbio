function theMonochromaticOptics = monochromaticVlambdaWeightedOpticsFromOI(theOI)

    theOptics = oiGet(theOI,'optics');
    thePSF = opticsGet(theOptics,'psf data');

    psfSupport = opticsGet(theOptics,'psf support','cyclesperdeg');
    spatialSupportXGridArcMin = psfSupport{1}*60;
    spatialSupportArcMin = spatialSupportXGridArcMin(1,:);
    wavelengthSupport = opticsGet(theOptics,'wave');

    micronsPerDegree = oiGet(theOI,'optics dist per deg')*1e6;
    pupilDiameterMM = opticsGet(theOptics,'pupildiameter')*1e3;

    % Load vLambda
    load T_xyz1931;
    vLambda = SplineCmf(S_xyz1931, T_xyz1931(2,:), WlsToS(wavelengthSupport));
    vLambda = vLambda / sum(vLambda);

    % Weight PSF by vLambda
    for iW = 1:numel(vLambda)
        if (iW == 1)
            theVLambdaWeightedPSF = squeeze(thePSF(:,:,iW)) * vLambda(iW);
        else
            theVLambdaWeightedPSF = theVLambdaWeightedPSF + squeeze(thePSF(:,:,iW)) * vLambda(iW);
        end
    end

    theMonochromaticPSF = thePSF * 0;
    for iW = 1:numel(vLambda)
        theMonochromaticPSF(:,:,iW) = theVLambdaWeightedPSF;
    end

    % Generate optics that has the same monochromatic PSF at each and every wavelength
    theMonochromaticOptics = oiFromPSF(theMonochromaticPSF, wavelengthSupport, ...
        spatialSupportArcMin, pupilDiameterMM, micronsPerDegree);
end
