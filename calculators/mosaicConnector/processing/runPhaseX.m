function runPhaseX(runParams)

    sfsExamined = [0.5 1.0 2.0 4.0 8.0 16.0];
    testSpatialFrequencyCPD = 4.0;
    photocurrentResponseDataFile = sprintf('LconeIsolatingResponses_%2.1fCPD', testSpatialFrequencyCPD);
    
    recomputePhotocurrents = ~true;
    if (recomputePhotocurrents)
        testLMScontrast = [0.1 0.0 0.0];
        instancesNum = 2;
        stimDurationSeconds = 0.5;
        fprintf('Will compute %d instances, each %2.1f seconds long\n', ...
            instancesNum, stimDurationSeconds);
        computePhotocurrents(runParams, instancesNum, stimDurationSeconds, ...
            testSpatialFrequencyCPD, testLMScontrast, photocurrentResponseDataFile);
    else
        computeRGCresponses(runParams, photocurrentResponseDataFile);
    end
end

function computeRGCresponses(runParams, photocurrentResponseDataFile)

    mFile = matfile(sprintf('%s.mat',photocurrentResponseDataFile), 'Writable', false);
    
    fprintf('\nImporting the cone mosaic ...');
    theConeMosaic = mFile.theConeMosaic;
    fprintf('Done !\n');
    
    fprintf('\nImporting the RGC mosaic ...');
    theMidgetRGCmosaic = mFile.theMidgetRGCmosaic;
    fprintf('Done !\n');
    
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'lineMarkerSize', 6, ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 30, ...
            'figureHeightInches', 14);
        
    % Visualize the mosaics
    visualizeMosaics = true;
    if (visualizeMosaics)
        fprintf('\nVisualizing the RGC mosaic with the optical image ...');
        theOISequence = mFile.theOIsequence;
        theFirstOI = theOISequence.frameAtIndex(1);
        zLevels = [0.3 1];
        hFig = visualizeConeAndRGCmosaicsWithRetinalImage(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
            theMidgetRGCmosaic, zLevels, 'centers', theFirstOI);
        fprintf('Done !\n');
    end
    plotlabOBJ.exportFig(hFig, 'pdf', sprintf('%s.mat',photocurrentResponseDataFile), pwd());


    thePresynapticResponses = mFile.isomerizations;
    %thePresynapticResponses = mFile.photocurrents;
    
    % Compute responses
    fprintf('\nComputing RGC responses ...');
    centerResponses = computeSubregionResponses(theMidgetRGCmosaic.centerWeights, thePresynapticResponses);
    surroundResponses = computeSubregionResponses(theMidgetRGCmosaic.surroundWeights, thePresynapticResponses);
    fprintf('Done !\n');
    
    % Normalize separately for center/surround
    centerResponses = centerResponses / max(centerResponses(:));
    surroundResponses = surroundResponses / max(surroundResponses(:));
    
    % Compute time axis
    timeAxis = (1:size(centerResponses,3))*theConeMosaic.integrationTime;
    
%     hFig = figure(98); clf;
%     plotWidthMicrons = 30;
%     plotHeightMicrons = 25;
%     
%     rgcsNum = size(centerResponses,2);
%     for iRGC = 1:rgcsNum
%         ax = axes('Position', [rgcRFpos(iRGC,1) rgcRFpos(iRGC,2) plotWidthMicrons plotHeightMicrons]/);
%         line(ax, timeAxis, squeeze(centerResponses(:,iRGC,:)), 'Color', [1 0 0]); hold(ax, 'on');
%         line(ax, timeAxis, squeeze(surroundResponses(:,iRGC,:)), 'Color', [0 0 1]);
%         set(ax, 'XTick', (0:50:500)/1000, 'YLim', [0 max(centerResponses(:))], 'YTick', 0:0.2:1);
%         set(ax, 'XTickLabel', {}, 'YTickLabel', {}, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');
%     end
end

function responses = computeSubregionResponses(weights, presynapticResponses)

    % Get dimensionalities
    [instancesNum, conesNum, timeBins] = size(presynapticResponses);
    rgcsNum = size(weights,2);
    
    % Form response matrix
    responses = zeros(instancesNum, rgcsNum, timeBins);
    for instanceIndex = 1:instancesNum
        % All presynaptic spatiotemporal responses for this instance
        instancePresynapticResponse = squeeze(presynapticResponses(instanceIndex,:,:));
        parfor iRGC = 1:rgcsNum
            % The RGC's weights
            iRGCweights = full(squeeze(weights(:,iRGC)));
            % The RGC temporal response
            responses(instanceIndex,iRGC,:) = iRGCweights' * instancePresynapticResponse;
        end
    end
    % mean over instances
    responses = mean(responses,1);
end



function computePhotocurrents(runParams, instancesNum, stimDurationSeconds, testSpatialFrequencyCPD, testLMScontrast, photocurrentResponseDataFile)
    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
    
    % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
    mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
    
    % Generate cone mosaic and connected mRGC mosaic patches
    [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams);
    
    % Generate optics appropriate for the RGC mosaic eccentricity
    wavelengthSampling = theConeMosaic.pigment.wave;
    pupilDiameterMM = 3.0;
    wavelengthsListToCompute = wavelengthSampling;
    wavefrontSpatialSamples = 1001;
    micronsPerDegree = []; % empty so as to compute for each eccentricity
    imposedRefractionErrorDiopters = 0;
    deltaEcc = 1;
    eccXrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(1))*[1 1];
    eccYrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(2))*[1 1];
    PolansSubjectID = 9;
    
    [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, wavefrontSpatialSamples, ...
                eccXrangeDegs, eccYrangeDegs, deltaEcc);

    theOI = theOIs{1,1,1};
    visualizePSFs(theOI, eccXrangeDegs(1), eccYrangeDegs(1));
    
    % Generate sine-wave scene as realized on a display
    stimColor = struct(...
        'backgroundChroma', [0.31, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', testLMScontrast);
    
    stimTemporalParams = struct(...
        'temporalFrequencyHz', 4.0, ...
        'stimDurationSeconds', stimDurationSeconds);
    
    stimSpatialParams = struct(...
        'fovDegs', max(theConeMosaic.fov),...
        'pixelsNum', 256, ...
        'gaborPosDegs', [0 0], ...
        'gaborSpatialFrequencyCPD', testSpatialFrequencyCPD, ...
        'gaborSigma', Inf, ...
        'gaborOrientationDegs', 0, ...
        'deltaPhaseDegs', 30);
    
    % Generate scenes corresponding to each spatial phase of a drifting grating
    theSceneFrames = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling);
    
    % Generate the OISequence
    fprintf('\nComputing the oiSequence ...');
    tic
    theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOI, stimSpatialParams, stimTemporalParams);
    %theOIsequence.visualize('montage', 'backendrenderer', 'figure');
    fprintf('Done in %2.1f minutes!\n', toc/60);
    
    % Generate eye movements for the oiSequence
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    emPaths = zeros(instancesNum, eyeMovementsNum, 2);
        
    % Compute isomerizations and photocurrents
    fprintf('\nComputing the mosaic response ...');
    tic
    
    [isomerizations, photocurrents, osLinearFilters] = ...
        theConeMosaic.computeForOISequence(theOIsequence, ...
        'emPaths', emPaths, 'currentFlag', true, ...
        'workerID', 1);
    
    fprintf('Done in %2.1f minutes!\n', toc/60);
    
    fprintf('\nExporting data ...');
    tic
    save(sprintf('%s.mat',photocurrentResponseDataFile), ...
        'theConeMosaic', ...
        'theMidgetRGCmosaic', ...
        'theOIsequence', ...
        'isomerizations', ...
        'photocurrents', ...
        '-v7.3');
    fprintf('Done in %2.1f minutes!\n', toc/60);
end

