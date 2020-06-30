function runPhaseX(runParams)

    % Location where to save all response files
    saveDir = strrep(fileparts(which(mfilename())), 'processing', 'responseFiles');
    
    % Whether to recompute the mosaics. This has to be done when we switch
    % to a new eccentricity
    recomputeMosaicsAndOptics = ~true;
    [theConeMosaic, theMidgetRGCmosaic, theOI] = mosaicsAndOpticsForEccentricity(runParams, recomputeMosaicsAndOptics, saveDir);
    
    testSpatialFrequenciesCPD = logspace(log10(0.1), log10(10),10);
    testTemporalFrequency = 2.0;
    instancesNum = 10;
    stimDurationSeconds = 1.0;
        
    
    testLMScontrast = [0.5 0.5 0.5];
    coneResponsesFileName = 'Achromatic50percentResponses';
    
    recomputePhotocurrents = true;
    recomputeNullResponses = true;
    
    if (recomputePhotocurrents)
        computeConeResponses(runParams, theConeMosaic, theOI, ...
            recomputeNullResponses, ...
            instancesNum, stimDurationSeconds, ...
            testSpatialFrequenciesCPD, testTemporalFrequency, testLMScontrast, ...
            coneResponsesFileName, saveDir);
    else
        computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
            testSpatialFrequenciesCPD, coneResponsesFileName, saveDir);
    end
end

function computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
    testSpatialFrequenciesCPD, coneResponsesFileName, saveDir)

    presynapticSignal = 'isomerizations';
    %presynapticSignal = 'photocurrents';
      
    mFile = matfile(fullfile(saveDir,nullResponseFilename(runParams)), 'Writable', false);
    if (strcmp(presynapticSignal,  'isomerizations'))
        theNullPresynapticResponses = mFile.isomerizationsNull;
    else
        theNullPresynapticResponses = mFile.photocurrentsNull;
    end


    for sfIndex = 1:numel(testSpatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = testSpatialFrequenciesCPD(sfIndex);
        dataFile = fullfile(saveDir,sprintf('%s_%2.1fCPD.mat',coneResponsesFileName, gaborSpatialFrequencyCPD));
        mFile = matfile(dataFile, 'Writable', false);
    
        % The stimulus sequence
        % theStimulusRGBsequence = mFile.theStimulusRGBsequence
        
        if (strcmp(presynapticSignal,  'isomerizations'))
            thePresynapticResponses = mFile.isomerizations;
        else
            thePresynapticResponses = mFile.photocurrents;
        end

        % Generate differential response
        thePresynapticResponses = bsxfun(@minus, thePresynapticResponses, theNullPresynapticResponses);

        % Compute responses
        fprintf('\nComputing RGC responses ...');
        tic
        centerResponses(sfIndex,:,:,:) = computeSubregionResponses(theMidgetRGCmosaic.centerWeights, thePresynapticResponses);
        surroundResponses(sfIndex,:,:,:) = -computeSubregionResponses(theMidgetRGCmosaic.surroundWeights, thePresynapticResponses);
        fprintf('Done in %2.1f minutes\n', toc/60);

        % Compute time axis
        timeAxis = (1:size(centerResponses,4))*theConeMosaic.integrationTime;
    end
    
    integratedResponses = centerResponses + surroundResponses;
    
    % Average over iterations
    centerResponsesMean = squeeze(mean(centerResponses,2));
    surroundResponsesMean = squeeze(mean(surroundResponses, 2));
    integratedResponsesMean = squeeze(mean(integratedResponses,2));
    
    maxMeanCenterResponse = max(abs(centerResponsesMean(:))); 
    maxMeanSurroundResponse = max(abs(surroundResponsesMean(:)));
    maxMeanIntegratedResponse = max(abs(integratedResponsesMean(:)));
    
    % Normalize separately for center/surround
    centerResponsesMean = centerResponsesMean / maxMeanCenterResponse;
    surroundResponsesMean = surroundResponsesMean / maxMeanCenterResponse;
    integratedResponsesMean = integratedResponsesMean / maxMeanCenterResponse;
    
    % Compute positions of RGCs
    RGCpositions = determineResponsePositions(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, theMidgetRGCmosaic.centerWeights);
    
    rgcsNum = size(RGCpositions,1);
    w = 0.1;
    h = 0.1;
    gw = 1-w;
    gh = 1-h;
    m = 0.03;
    for sfIndex = 1:numel(testSpatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = testSpatialFrequenciesCPD(sfIndex);
        hFig = figure(100+sfIndex); clf;
        set(hFig, 'name', sprintf('%2.1f c/deg', gaborSpatialFrequencyCPD), 'Position', [10 10 1100 1000]);
        for iRGC = 1:rgcsNum
            ax = axes('Position', [gw*RGCpositions(iRGC,1)+m gh*RGCpositions(iRGC,2)+m w h]);
            line(ax, timeAxis, timeAxis*0, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
            hold(ax, 'on');
            %line(ax, timeAxis, squeeze(centerResponses(sfIndex,:,iRGC,:)), 'Color', [1 0 0], 'LineWidth', 1.0);
            %line(ax, timeAxis, squeeze(surroundResponses(sfIndex,:,iRGC,:)), 'Color', [0 0 1], 'LineWidth', 1.0);
            line(ax, timeAxis, squeeze(centerResponsesMean(sfIndex,iRGC,:)), 'Color', [1 0 0], 'LineWidth', 1.0);
            line(ax, timeAxis, squeeze(surroundResponsesMean(sfIndex,iRGC,:)), 'Color', [0 0 1], 'LineWidth', 1.0);
            line(ax, timeAxis, squeeze(integratedResponsesMean(sfIndex,iRGC,:)), 'Color', [0 0 0], 'LineWidth', 1.5);
            set(ax, 'XTickLabel', {}, 'YTickLabel', {}, 'XTick', (0:50:1000)/1000, 'YTick', -1:0.25:1, ...
                'XLim', [0 timeAxis(end)], 'YLim', [-1 1]);
            box(ax, 'on'); grid(ax, 'on');
        end
    end
    
%     plotlabOBJ = plotlab();
%     plotlabOBJ.applyRecipe(...
%             'colorOrder', [1 0 0; 0 0 1], ...
%             'axesBox', 'off', ...
%             'axesTickDir', 'in', ...
%             'renderer', 'painters', ...
%             'lineMarkerSize', 6, ...
%             'axesTickLength', [0.01 0.01], ...
%             'legendLocation', 'SouthWest', ...
%             'figureWidthInches', 30, ...
%             'figureHeightInches', 14);
        
    % Visualize the mosaics
    visualizeMosaics = ~true;
    if (visualizeMosaics)
        fprintf('\nVisualizing the RGC mosaic with the optical image ...');
        theOISequence = mFile.theOIsequence;
        theFirstOI = theOISequence.frameAtIndex(1);
        zLevels = [0.3 1];
        hFig = visualizeConeAndRGCmosaicsWithRetinalImage(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
            theMidgetRGCmosaic, zLevels, 'centers', theFirstOI); 
        plotlabOBJ.exportFig(hFig, 'pdf', sprintf('%s.mat',coneResponsesFileName), pwd());
        fprintf('Done !\n');
    end
    


   
    
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

function RGCpositions = determineResponsePositions(theConeMosaic, eccentricityMicrons, centerWeights)
     % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
    
    rgcsNum = size(centerWeights,2);
    RGCpositions = zeros(rgcsNum,2);
    for iRGC = 1:rgcsNum
        weights = full(squeeze(centerWeights(:, iRGC)));
        centerIndices = find(weights>0);
        RGCpositions(iRGC,:) = mean(conePositionsMicrons(centerIndices,:),1);
    end
    
    % Normalize to [0..1]
    mins = min(RGCpositions, [], 1);
    maxs = max(RGCpositions, [], 1);
    RGCpositions(:,1) = (RGCpositions(:,1)-mins(1))/(maxs(1)-mins(1));
    RGCpositions(:,2) = (RGCpositions(:,2)-mins(2))/(maxs(2)-mins(2));
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
end

function [theConeMosaic, theMidgetRGCmosaic, theOI] = mosaicsAndOpticsForEccentricity(runParams, recompute, saveDir)

    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
        
    mosaicsAndOpticsFileName = sprintf('MosaicsAndOpticsForEccentricity_%2.0f_%2.0f_%2.0f_%2.0f_microns_coneSpecificity_%2.0f_orphanPolicy_%s.mat', ...
        mosaicParams.rgcMosaicPatchEccMicrons(1), mosaicParams.rgcMosaicPatchEccMicrons(2), ...
        mosaicParams.rgcMosaicPatchSizeMicrons(1), mosaicParams.rgcMosaicPatchSizeMicrons(2), ...
        mosaicParams.maximizeConeSpecificity, mosaicParams.orphanRGCpolicy);
    
    tic  
    if (recompute)
        fprintf('\nComputing mosaics and optics ...');
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
        save(fullfile(saveDir,mosaicsAndOpticsFileName), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOI');
    else
        fprintf('\nLoading mosaics and optics ...');
        load(fullfile(saveDir,mosaicsAndOpticsFileName), 'theConeMosaic', 'theMidgetRGCmosaic', 'theOI')
    end
    fprintf('Done in %2.1f minutes\n', toc/60);
    
    displayPSFs = false;
    if (displayPSFs)
        visualizePSFs(theOI, eccXrangeDegs(1), eccYrangeDegs(1));
    end
end

function computeConeResponses(runParams, theConeMosaic, theOI, ...
    recomputeNullResponses, ...
    instancesNum, stimDurationSeconds, ...
    testSpatialFrequenciesCPD, testTemporalFrequency, testLMScontrast, ...
    coneResponsesFileName, saveDir)

    % Generate sine-wave scene as realized on a display
    stimColor = struct(...
        'backgroundChroma', [0.31, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', testLMScontrast);
    
    stimTemporalParams = struct(...
        'temporalFrequencyHz',testTemporalFrequency, ...
        'stimDurationSeconds', stimDurationSeconds);
    
    stimSpatialParams = struct(...
        'fovDegs', max(theConeMosaic.fov),...
        'pixelsNum', 256, ...
        'gaborPosDegs', [0 0], ...
        'gaborSpatialFrequencyCPD', 0, ...
        'gaborSigma', Inf, ...
        'gaborOrientationDegs', 0, ...
        'deltaPhaseDegs', 30);
    
    if (recomputeNullResponses)
        % Compute the null (zero contrast, no-noise responses)
        nullColor = stimColor;
        nullColor.lmsContrast = [0.0 0.0 0.0];
        wavelengthSampling = theConeMosaic.pigment.wave;
        theNullSceneFrames = generateStimulusFrames(nullColor, stimSpatialParams, wavelengthSampling);
        fprintf('\nComputing the NULL oiSequence ...');
        tic
        theOIsequence = generateOISequenceForDriftingGrating(theNullSceneFrames, theOI, stimSpatialParams, stimTemporalParams);
        fprintf('Done in %2.1f minutes\n', toc/60);
        eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
        emPaths = zeros(1, eyeMovementsNum, 2);

        % Save original isomerizations and photocurrent noise flags
        originalIsomerizationNoiseFlag = theConeMosaic.noiseFlag;
        originalPhotocurrentNoiseFlag = theConeMosaic.os.noiseFlag;

        % Set noiseFlags to none for the null response
        theConeMosaic.noiseFlag = 'none';
        theConeMosaic.os.noiseFlag = 'none';

        fprintf('\nComputing the NULL mosaic responses ...');
        tic
        % Compute null responses, photocurrent impulse responses and mean photocurrents
        [isomerizationsNull, photocurrentsNull, osImpulseResponseFunctions, osMeanCurrents] = ...
                theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', emPaths, 'currentFlag', true);
        fprintf('Done in %2.1f minutes\n', toc/60);
        
        % Restore noise flags
        theConeMosaic.noiseFlag = originalIsomerizationNoiseFlag;
        theConeMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;

        for k = 1:theOIsequence.length
            theFrame = xyz2rgb(oiGet(theOIsequence.frameAtIndex(k), 'xyz'));
            if (k == 1)
                theStimulusRGBsequence = zeros(theOIsequence.length, size(theFrame,1), size(theFrame,2), size(theFrame,3));
            end
            theStimulusRGBsequence(k,:,:,:) = theFrame;
        end
        
        fprintf('\nSaving the NULL responses ...');
        tic
        % Save the null responses, and the photocurrent impulse responses/mean currents
        save(fullfile(saveDir,nullResponseFilename(runParams)), ...
             'theStimulusRGBsequence', ...
             'isomerizationsNull', ...
             'photocurrentsNull', ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents', ...
             '-v7.3');
         fprintf('Done in %2.1f minutes\n', toc/60);
    else
        fprintf('\nLoading null responses ...');
        tic
        % Load the photocurrent impulse responses/mean currents
        load(fullfile(saveDir,nullResponseFilename(runParams)), ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents');
         fprintf('Done in %2.1f minutes\n', toc/60);
    end
    
    % Compute test responses
    for sfIndex = 1:numel(testSpatialFrequenciesCPD)
        % Test spatial frequency
        stimSpatialParams.gaborSpatialFrequencyCPD = testSpatialFrequenciesCPD(sfIndex);
        
        % Generate scenes corresponding to each spatial phase of a drifting grating
        wavelengthSampling = theConeMosaic.pigment.wave;
        theSceneFrames = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling);

        % Generate the OISequence
        fprintf('\nComputing the %2.1f c/deg, [%2.1f %2.1f %2.1f] LMS %2.1f second stimulus oiSequence ...', ...
            stimSpatialParams.gaborSpatialFrequencyCPD, ...
            testLMScontrast(1), testLMScontrast(2), testLMScontrast(3), ...
            stimTemporalParams.stimDurationSeconds);
        tic
        theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOI, stimSpatialParams, stimTemporalParams);
        %theOIsequence.visualize('montage', 'backendrenderer', 'figure');
        fprintf('Done in %2.1f minutes!\n', toc/60);
    
        % Generate eye movements for the oiSequence
        eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
        emPaths = zeros(instancesNum, eyeMovementsNum, 2);

        % Compute isomerizations and photocurrents
        fprintf('\nComputing %d response instances ...', instancesNum);
        tic
        [isomerizations, photocurrents] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...
            'interpFilters', osImpulseResponseFunctions, ...   % employ photocurrent impulse responses from null response
            'meanCur', osMeanCurrents, ...                     % employ mean photocurrents from null response
            'currentFlag', true);
    
        fprintf('Done in %2.1f minutes!\n', toc/60);
    
        % Save data
        fprintf('\nExporting the %2.1f c/deg mosaic responses ...', stimSpatialParams.gaborSpatialFrequencyCPD);
        tic
        for k = 1:theOIsequence.length
            theFrame = xyz2rgb(oiGet(theOIsequence.frameAtIndex(k), 'xyz'));
            if (k == 1)
                theStimulusRGBsequence = zeros(theOIsequence.length, size(theFrame,1), size(theFrame,2), size(theFrame,3));
            end
            theStimulusRGBsequence(k,:,:,:) = theFrame;
        end
        save(fullfile(saveDir,sprintf('%s_%2.1fCPD.mat',coneResponsesFileName, stimSpatialParams.gaborSpatialFrequencyCPD)), ...
            'theStimulusRGBsequence', ...
            'isomerizations', ...
            'photocurrents', ...
            '-v7.3');
        fprintf('Done in %2.1f minutes!\n', toc/60);
    end % sfIndex
    
end

function fname = nullResponseFilename(runParams)
    fname = sprintf('NullResponses_%2.0f_%2.0f_%2.0f_%2.0f.mat', ...
            runParams.rgcMosaicPatchEccMicrons(1), runParams.rgcMosaicPatchEccMicrons(2), ...
            runParams.rgcMosaicPatchSizeMicrons(1), runParams.rgcMosaicPatchSizeMicrons(2));
end
        