function [theConeMosaic, theMidgetRGCmosaic, theOptics] = mosaicsAndOpticsForEccentricity(runParams, recompute, saveDir)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
        
    tic  
    if (recompute)
        fprintf('\nComputing mosaics and optics ...');
        % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
        mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));

        % Generate cone mosaic and connected mRGC mosaic patches
        [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams, ...
             runParams.outputFile, runParams.exportsDir);
         
        wavelengthSampling = theConeMosaic.pigment.wave;
        pupilDiameterMM = 3.0;
        wavelengthsListToCompute = wavelengthSampling;
        wavefrontSpatialSamples = 501;
        micronsPerDegree = []; % empty so as to compute for each eccentricity
        imposedRefractionErrorDiopters = 0;
        deltaEcc = 1;
        eccXrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(1))*[1 1];
        eccYrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(2))*[1 1];

        [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(runParams.PolansWavefrontAberrationSubjectID, ...
                    imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, ...
                    wavefrontSpatialSamples, eccXrangeDegs, eccYrangeDegs, deltaEcc);

        theOptics = theOIs{1,1,1};
        save(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOptics');
    else
        fprintf('\nLoading mosaics and optics ...');
        load(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOptics')
    end

    fprintf('Done in %2.1f minutes\n', toc/60);
    
    displayPSFs = false;
    if (displayPSFs)
        visualizePSFs(theOI, eccXrangeDegs(1), eccYrangeDegs(1));
    end
end
