function [theConeMosaic, theMidgetRGCmosaic, theOptics] = mosaicsAndOpticsForEccentricity(runParams, recomputeMosaicAndOptics, recomputeOpticsOnly, saveDir)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
        
    % Extract the Polans subjectID for the optics
    runParams.PolansSubjectID = runParams.deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage;
    if (numel(runParams.PolansSubjectID )>1)
        runParams.PolansSubjectID = runParams.PolansSubjectID(1);
        fprintf(2, 'deconvolutionOpticsParams indicates more than 1 subject being used to derive the deconvolution model. Will generate optics for the first of these subjects\n');
    end
        
    tic  
    if (recomputeMosaicAndOptics)
        fprintf('\nComputing mosaics and optics ...');
        % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
        mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
       
        % Generate cone mosaic and connected mRGC mosaic patches
        [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams, ...
             runParams.deconvolutionOpticsParams, runParams.outputFile, runParams.exportsDir);
        
        % Compute the optics
        [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(runParams.PolansSubjectID, runParams.noLCA, ...
            theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);
        
        % Save mosaic and optics
        save(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOptics');
        
    elseif (recomputeOpticsOnly)
        % Load previously generated mosaic
        fprintf('\nLoading previously-generated mosaics ...');
        load(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic');
        
        fprintf('\nComputing new optics with noLCA flag: %d\n', runParams.noLCA);
        [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(runParams.PolansSubjectID, runParams.noLCA, ...
            theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);
        
        % Save mosaic and optics
        save(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOptics');
    else
        fprintf('\nLoading previously-generated mosaics & optics ...');
        load(fullfile(saveDir,mosaicsAndOpticsFileName(runParams)), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOptics')
    end

    fprintf('Done in %2.1f minutes\n', toc/60);
    
    
end

function [theOptics, eccXrangeDegs, eccYrangeDegs] =  generatePolansOptics(PolansSubjectID, noLCA, wavelengthSampling, rgcMosaicPatchEccMicrons)
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
                'noLCA', noLCA);

    theOptics = theOIs{1,1,1};
end
