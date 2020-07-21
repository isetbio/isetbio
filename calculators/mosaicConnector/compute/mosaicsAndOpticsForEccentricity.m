function [theConeMosaic, theMidgetRGCmosaic, theOptics, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsForEccentricity(runParams, recomputeMosaicAndOptics, recomputeOpticsOnly, saveDir)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
        
    tic  
    
    % Compute filenames for mosaics and optics
    [mosaicsFilename, opticsFilename, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
    
    if (recomputeMosaicAndOptics)
        fprintf('\nComputing mosaics and optics ...');
        % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
        mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
       
        % Generate cone mosaic and connected mRGC mosaic patches
        [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams, ...
             runParams.deconvolutionOpticsParams, runParams.outputFile, runParams.exportsDir);
        
        % Compute the optics
        fprintf('\nComputing optics with noLCA flag: %d and noOptics flag: %d\n', runParams.noLCA, runParams.noOptics);
        [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(PolansSubjectID, runParams.noLCA, runParams.noOptics, ...
            theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);
        
        % Save mosaic and optics
        save(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theMidgetRGCmosaic');
        save(fullfile(saveDir,opticsFilename),'theOptics');
        
    elseif (recomputeOpticsOnly)
        % Load previously generated mosaic
        fprintf('\nLoading previously-generated mosaics ...');
        load(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theMidgetRGCmosaic');
        
        fprintf('\nComputing new optics with noLCA flag: %d and noOptics flag: %d\n', runParams.noLCA, runParams.noOptics);
        [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(PolansSubjectID, runParams.noLCA, runParams.noOptics, ...
            theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);
        
        % Save mosaic and optics
        save(fullfile(saveDir,opticsFilename),'theOptics');
        
    else
        fprintf('\nLoading previously-generated mosaics & optics ...');
        load(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theMidgetRGCmosaic');
        load(fullfile(saveDir,opticsFilename),'theOptics');
    end

    fprintf('Done in %2.1f minutes\n', toc/60);
    
    
end

function [theOptics, eccXrangeDegs, eccYrangeDegs] =  generatePolansOptics(PolansSubjectID, noLCA, noOptics, wavelengthSampling, rgcMosaicPatchEccMicrons)
    pupilDiameterMM = 3.0;
    wavelengthsListToCompute = wavelengthSampling;
    wavefrontSpatialSamples = 501;
    micronsPerDegree = []; % empty so as to compute for each eccentricity
    imposedRefractionErrorDiopters = 0;
    deltaEcc = 1;
    eccXrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*rgcMosaicPatchEccMicrons(1))*[1 1];
    eccYrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*rgcMosaicPatchEccMicrons(2))*[1 1];

    [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, ...
                wavefrontSpatialSamples, eccXrangeDegs, eccYrangeDegs, deltaEcc, ...
                'noLCA', noLCA, 'noOptics', noOptics);

    theOptics = theOIs{1,1,1};
end
