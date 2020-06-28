function runPhaseX(runParams)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
    
    % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
    mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
    
    % Generate cone mosaic and connected mRGC mosaic patches
    %[theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams);
    
    % Generate optics appropriate for the RGC mosaic eccentricity
    pupilDiameterMM = 3.0;
    wavelengthsListToCompute = [400:50:750];
    wavefrontSpatialSamples = 1001;
    micronsPerDegree = []; % empty so as to compute for each eccentricity
    imposedRefractionErrorDiopters = 0;
    deltaEcc = 1;
    eccXrangeDegs = 0*WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(1))*[1 1];
    eccYrangeDegs = 0*WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(2))*[1 1];
    PolansSubjectID = 9;
    [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, wavefrontSpatialSamples, ...
                eccXrangeDegs, eccYrangeDegs, deltaEcc);
    theOI = theOIs{1,1,1};
    visualizePSFs(theOI, eccXrangeDegs(1), eccYrangeDegs(1));
    
    % Generate sine-wave scene
    
    % Compute static cone mosaic isomerizations
    
    % Compute static mRGC RF responses
    
end

function visualizePSFs(theOI, eccX, eccY)
    optics = oiGet(theOI, 'optics');
    wavelengthSupport = opticsGet(optics, 'wave');
    thePSFs = opticsGet(optics,'psf data');
    
    % Extract support in arcmin
    psfSupportMicrons = opticsGet(optics,'psf support','um');
    if (isfield(optics, 'micronsPerDegree'))
        micronsPerDegree = optics.micronsPerDegree;
    else
        focalLengthMeters = opticsGet(optics, 'focalLength');
        focalLengthMicrons = focalLengthMeters * 1e6;
        micronsPerDegree = focalLengthMicrons * tand(1);
    end

    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    thePSFsupportDegs = xGridDegs(1,:);
    
    xRange = 5/60*[-1 1];
    figure(1); clf;
    for wIndex = 1:numel(wavelengthSupport)
        subplot(2,4,wIndex);
        imagesc(thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(:,:,wIndex)));
        axis 'square'
        set(gca, 'XLim', xRange, 'YLim', xRange);
        title(sprintf('%2.0f nm (%2.1f, %2.1f) degs', wavelengthSupport(wIndex), eccX, eccY));
    end
    colormap(gray);
end

