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
    wavelengthsListToCompute = [550];
    wavefrontSpatialSamples = 1001;
    micronsPerDegree = []; % empty so as to compute for each eccentricity
    imposedRefractionErrorDiopters = 0;
    deltaEcc = 1;
    eccXrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(1))*[1 1];
    eccYrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(2))*[1 1];
    PolansSubjectID = 8;
    [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, wavefrontSpatialSamples, ...
                eccXrangeDegs, eccYrangeDegs, deltaEcc);
    theOI = theOIs{1,1,1};
    
    
    % Generate sine-wave scene
    
    % Compute static cone mosaic isomerizations
    
    % Compute static mRGC RF responses
    
end


