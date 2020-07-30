function [theConeMosaic, theMidgetRGCmosaic, theOptics, opticsPostFix, PolansSubjectID] = ...
    mosaicsAndOpticsForEccentricity(runParams, recomputeConeMosaic, recomputeRGCmosaic, recomputeOptics, saveDir, ...
    visualizeCronerKaplanDeconvolutionModel)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
        
    tic  
    
    % Compute filenames for mosaics and optics
    [mosaicsFilename, opticsFilename, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
    
    % Determine what we need to recompute and what can be loaded from the existing files
    recomputeOpticsOnly = ((~recomputeConeMosaic) && (~recomputeRGCmosaic)) && (recomputeOptics);
    
    if (recomputeConeMosaic) || (recomputeRGCmosaic)
        
        if (~recomputeConeMosaic)
            % Load previously generated mosaic
            fprintf('\nLoading previously-generated cone mosaic ...');
            load(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theConeMosaicMetaData');
        else
            fprintf('\nRecomputing cone mosaic ...');
            % Generate a regular hex cone mosaic patch with the desired eccentricity and size
            extraMicronsForSurroundCones = estimateMaxSurroundRadiusMicrons(mosaicParams.rgcMosaicPatchEccMicrons, mosaicParams.rgcMosaicPatchSizeMicrons, runParams.deconvolutionOpticsParams);
            [theConeMosaic, coneMosaicEccDegs, coneMosaicSizeMicrons, conePositionsMicrons, coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones] = ...
                generateRegularHexMosaicPatch(...
                mosaicParams.rgcMosaicPatchEccMicrons, ...
                mosaicParams.rgcMosaicPatchSizeMicrons, ...
                extraMicronsForSurroundCones);
            
            % Assemble metadata struct
            theConeMosaicMetaData = struct(...
                'coneMosaicEccDegs', coneMosaicEccDegs, ...
                'coneMosaicSizeMicrons', coneMosaicSizeMicrons, ...
                'conePositionsMicrons', conePositionsMicrons, ...
                'coneSpacingsMicrons', coneSpacingsMicrons, ...
                'coneTypes', coneTypes, ...
                'extraMicronsForSurroundCones', extraMicronsForSurroundCones ...
                );
        end
        
        fprintf('\nRecomputing mRGC mosaic and connecting it to the cone mosaic ...');
        % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
        mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
       
        % Generate connected mRGC mosaic patch
        theMidgetRGCmosaic = generateMRGCMosaicConnectedToConeMosaic(theConeMosaicMetaData, mRGCmosaicFile, mosaicParams, ...
             runParams.deconvolutionOpticsParams, visualizeCronerKaplanDeconvolutionModel, runParams.outputFile, runParams.exportsDir);
        
        % Save the mosaics
        save(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theConeMosaicMetaData', 'theMidgetRGCmosaic');
        
        if (recomputeOptics)
            % Generate the optics
            fprintf('\nComputing optics with noLCA flag: %d and noOptics flag: %d\n', runParams.noLCA, runParams.noOptics);
            [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(PolansSubjectID, runParams.noLCA, runParams.noOptics, ...
                theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);

            % Save the optics
            save(fullfile(saveDir,opticsFilename),'theOptics');
        else
            load(fullfile(saveDir,opticsFilename),'theOptics');
        end
        
    elseif (recomputeOpticsOnly)
        % Load previously generated mosaic
        fprintf('\nLoading previously-generated cone/mRGC mosaics ...');
        load(fullfile(saveDir,mosaicsFilename), 'theConeMosaic', 'theMidgetRGCmosaic');
        
        fprintf('\nComputing optics with noLCA flag: %d and noOptics flag: %d\n', runParams.noLCA, runParams.noOptics);
        [theOptics, eccXrangeDegs, eccYrangeDegs] = generatePolansOptics(PolansSubjectID, runParams.noLCA, runParams.noOptics, ...
            theConeMosaic.pigment.wave, mosaicParams.rgcMosaicPatchEccMicrons);
        
        % Save mosaic and optics
        save(fullfile(saveDir,opticsFilename),'theOptics');
        
    else
        fprintf('\nLoading previously-generated cone/mRGC mosaics & optics ...');
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

function [theConeMosaic, coneMosaicEccDegs, coneMosaicSizeMicrons, conePositionsMicrons, coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones] = ...
    generateRegularHexMosaicPatch(eccentricityMicrons, sizeMicrons, extraMicronsForSurroundCones)

    % Compute the cone mosaic FOV in degrees
    coneMosaicCenterPositionMM = eccentricityMicrons * 1e-3;
    coneMosaicSizeMicrons = sizeMicrons + 2*extraMicronsForSurroundCones*[1 1];
    coneMosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMM);
    fovDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(coneMosaicSizeMicrons, sqrt(sum((coneMosaicCenterPositionMM*1e3).^2,2.0)));
    
    % Determine the median cone spacing with the patch
    whichEye = 'right';
    coneSpacingMicrons = medianConeSpacingInPatch(whichEye, eccentricityMicrons, coneMosaicSizeMicrons);
    
    % Generate reg hex cone mosaic
    resamplingFactor = 5;
    spatialDensity = [0 0.5 0.25 0.15];
    sConeFreeRadiusMicrons = 0;
    theConeMosaic = coneMosaicHex(resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'micronsPerDegree', coneMosaicSizeMicrons(1)/fovDegs(1), ...
        'integrationTime', 5/1000, ...
        'customLambda', coneSpacingMicrons, ...
        'customInnerSegmentDiameter', coneSpacingMicrons * 0.7, ...
        'spatialDensity', spatialDensity, ...
        'sConeMinDistanceFactor', 2, ...
        'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );


    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
end

function extraMicronsForSurroundCones = estimateMaxSurroundRadiusMicrons(eccentricityMicrons, sizeMicrons, deconvolutionOpticsParams)
    w = WatsonRGCModel('generateAllFigures', false);
    posUnits = 'mm'; densityUnits = 'mm^2';
    coneEccMaxMicrons = max(abs([eccentricityMicrons+sizeMicrons/2; eccentricityMicrons-sizeMicrons/2]));
    
    % Cone spacing at the most eccentric position
    coneSpacingMicronsMax = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
        coneEccMaxMicrons*1e-3, 'right', posUnits, densityUnits, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
 
    % Compute RF params using the CronerKaplan model
    ck = CronerKaplanRGCModel('generateAllFigures', false, 'instantiatePlotLab', false);
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(coneSpacingMicronsMax, coneEccMaxMicrons, deconvolutionOpticsParams);
    retinalSurroundRadiusDegsAt1overE = synthesizedRFParams.retinal.surroundRadiiDegs;
    extraDegsForSurroundCones = retinalSurroundRadiusDegsAt1overE*2;
    extraMicronsForSurroundCones = ceil(WatsonRGCModel.sizeDegsToSizeRetinalMicrons(extraDegsForSurroundCones, synthesizedRFParams.eccDegs));
end