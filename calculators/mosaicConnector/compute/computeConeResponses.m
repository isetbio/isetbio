function computeConeResponses(runParams, ...
    stimColor,  stimTemporalParams, stimSpatialParams, ...
    theConeMosaic, theOptics, ...
    recomputeNullResponses, ...
    instancesNum, ...
    spatialFrequenciesCPD, ...
    saveDir)

    
    % Determine phase increment at each stimulus frame so as to get required temporal frequency with the
    % constraint that the frame duration is evenly divisible by cone mosaic integration time
    [stimSpatialParams.deltaPhaseDegs, stimTemporalParams.temporalFrequencyHz] = ...
        computeDeltaPhaseDegs(stimTemporalParams.temporalFrequencyHz, theConeMosaic.integrationTime);
    
    if (recomputeNullResponses)
        % Compute the null (zero contrast, no-noise responses)
        nullColor = stimColor;
        nullColor.lmsContrast = [0.0 0.0 0.0];
        
        % Generate the null stimulus frame scene sequence
        wavelengthSampling = theConeMosaic.pigment.wave;
        [theNullSceneFrames, presentationDisplay] = generateStimulusFrames(nullColor, stimSpatialParams, wavelengthSampling);
        
        % Generate the stimulus SRGBsequence
        [theStimulusSRGBsequence, ~, theStimulusLMSexcitationsSequence, ...
            stimulusSpatialSupportDegs, stimulusSpatialSupportMicrons] = extractSRGBsequence(theNullSceneFrames, presentationDisplay, theConeMosaic.qe);
        
        % Generate the corresponding optical image sequence
        fprintf('\nComputing the NULL oiSequence ...');
        tic

        theOIsequence = generateOISequenceForDriftingGrating(theNullSceneFrames, theOptics, ...
            stimTemporalParams);
        clear 'theNullSceneFrames';
        fprintf('Done in %2.1f minutes\n', toc/60);

        % Zero eye movements
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
        
        % plot the photocurrent impulse responses
        plotImpulseResponses = true;
        if (plotImpulseResponses)
            figure(222);
            plot(osImpulseResponseFunctions(:,1), 'r-', 'LineWidth', 1.5); hold on;
            plot(osImpulseResponseFunctions(:,2), 'g-', 'LineWidth', 1.5);
            plot(osImpulseResponseFunctions(:,3), 'b-', 'LineWidth', 1.5);
            title(sprintf('Mean currents: L=%2.2f, M=%2.2f, S=%2.2f', osMeanCurrents(1), osMeanCurrents(2), osMeanCurrents(3)));
        end
        
        % Restore noise flags
        theConeMosaic.noiseFlag = originalIsomerizationNoiseFlag;
        theConeMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;

        
        % Save the data
        fprintf('\nExporting the NULL responses ...');
        tic
        % Extract the RGB sequence of optical images and LMS excitation images
        [theOISRGBsequence, ~, theRetinalLMSexcitationsSequence, ...
            retinalSpatialSupportDegs, retinalSpatialSupportMicrons] = extractSRGBsequence(theOIsequence, presentationDisplay, theConeMosaic.qe);
        
        % Extract the stimulus and response time axes
        stimulusTimeAxis = theOIsequence.timeAxis;
        responseTimeAxis = theConeMosaic.timeAxis;
        
        % Save the null responses, and the photocurrent impulse responses/mean currents
        save(fullfile(saveDir, nullResponseFilename(runParams)), ...
             'stimColor', 'stimTemporalParams', 'stimSpatialParams', ...
             'theStimulusSRGBsequence', ...
             'theStimulusLMSexcitationsSequence', ...
             'theOISRGBsequence', ...
             'theRetinalLMSexcitationsSequence', ...
             'stimulusSpatialSupportDegs', 'stimulusSpatialSupportMicrons', ...
             'retinalSpatialSupportDegs', 'retinalSpatialSupportMicrons', ...
             'isomerizationsNull', ...
             'photocurrentsNull', ...
             'responseTimeAxis', ...
             'stimulusTimeAxis', ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents', ...
             '-v7.3');
         fprintf('Done in %2.1f minutes\n', toc/60);
         clear('theStimulusRGBsequence');
         clear('isomerizationsNull');
         clear('photocurrentsNull');
         clear('theOIsequence')
    else
        fprintf('\nLoading null responses ...');
        tic
        % Load the photocurrent impulse responses and mean currents
        load(fullfile(saveDir, nullResponseFilename(runParams)), ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents');
         fprintf('Done in %2.1f minutes\n', toc/60);
    end
    
    
    % Compute test responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        % Test spatial frequency
        stimSpatialParams.gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Generate the test (spatial frequency) stimulus frame scene sequence
        wavelengthSampling = theConeMosaic.pigment.wave;
        [theSceneFrames, presentationDisplay, sceneLuminanceSlice] = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling);
        
        % Generate the stimulus SRGBsequence
        [theStimulusSRGBsequence,~, theStimulusLMSexcitationsSequence, ...
            stimulusSpatialSupportDegs, stimulusSpatialSupportMicrons] = extractSRGBsequence(theSceneFrames, presentationDisplay, theConeMosaic.qe);
        
        figure(4000);
        subplot(3,5,sfIndex)
        imagesc(sceneLuminanceSlice);
        set(gca, 'CLim', [0 100]);
        maxSceneLuminance(sfIndex) = max(sceneLuminanceSlice(:));

        subplot(3,5,15);
        plot(spatialFrequenciesCPD(1:sfIndex), maxSceneLuminance(1:sfIndex), 'ks-');
        ylabel('scene luminance');
        xlabel('c/deg');
        set(gca, 'XLim', [0 100], 'XScale', 'log');
        colormap(gray);
        
        % Generate the corresponding optical image sequence
        fprintf('\nComputing the %2.1f c/deg, [%2.2f %2.2f %2.2f] LMS %2.1f-second stimulus oiSequence ...', ...
            stimSpatialParams.gaborSpatialFrequencyCPD, ...
            stimColor.lmsContrast(1), stimColor.lmsContrast(2), stimColor.lmsContrast(3), ...
            stimTemporalParams.stimDurationSeconds);
        tic
        theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOptics, stimTemporalParams);
        clear 'theSceneFrames';
        %theOIsequence.visualize('montage', 'backendrenderer', 'figure');
        fprintf('Done in %2.1f minutes!\n', toc/60);
        
        % Zero eye movements
        eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
        emPaths = zeros(instancesNum, eyeMovementsNum, 2);

        % Compute isomerizations and photocurrents using the meanCurrents
        % and outersegment impulse response functions computed from the
        % null stimulus
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
        fprintf('\nExporting the %2.1f c/deg responses ...', stimSpatialParams.gaborSpatialFrequencyCPD);
        tic
        % Extract the RGB sequence of optical images and LMS excitation images
        [theOISRGBsequence, retinalIlluminanceSlice, theRetinalLMSexcitationsSequence, ...
            retinalSpatialSupportDegs, retinalSpatialSupportMicrons] = extractSRGBsequence(theOIsequence, presentationDisplay, theConeMosaic.qe);
       
        
        figure(6000);
        subplot(3,5,sfIndex)
        imagesc(retinalIlluminanceSlice);
        set(gca, 'CLim', [1 2]);
        maxRetinalIlluminance(sfIndex) = max(retinalIlluminanceSlice(:));
        subplot(3,5,15);
        plot(spatialFrequenciesCPD(1:sfIndex), maxRetinalIlluminance(1:sfIndex), 'ks-');
        ylabel('retinal illuminance');
        xlabel('c/deg');
        set(gca, 'XLim', [0.1 60], 'XScale', 'log');
        colormap(gray);
        
        save(fullfile(saveDir,sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, stimColor.lmsContrast), stimSpatialParams.gaborSpatialFrequencyCPD)), ...
            'stimColor', 'stimTemporalParams', 'stimSpatialParams', ...
            'theStimulusSRGBsequence', ...
            'theStimulusLMSexcitationsSequence', ...
            'theOISRGBsequence', ...
            'theRetinalLMSexcitationsSequence', ...
            'stimulusSpatialSupportDegs', 'stimulusSpatialSupportMicrons', ...
            'retinalSpatialSupportDegs', 'retinalSpatialSupportMicrons', ...
            'isomerizations', ...
            'photocurrents', ...
            '-v7.3');
        
        fprintf('Done in %2.1f minutes!\n', toc/60);
        clear('theStimulusSRGBsequence');
        clear('isomerizations');
        clear('photocurrents');
        clear('theOIsequence');
        
    end % sfIndex
end

function [theSRGBsequence, theSlices, theLMSexcitationsSequence, spatialSupportDegs, spatialSupportMicrons] = extractSRGBsequence(theSequence, presentationDisplay, coneQuantalEfficiencies)
    displaySPDs = displayGet(presentationDisplay, 'spd'); 
    
    if (isa(theSequence, 'oiArbitrarySequence'))
        % Dealing with an optical image
        framesNum = theSequence.length;
        retinalIlluminanceImage = oiGet(theSequence.frameAtIndex(1), 'illuminance');
        midRow = round(size(retinalIlluminanceImage,1)/2);
         
        % Preallocate memory
        theSRGBsequence = zeros(framesNum , size(retinalIlluminanceImage,1), size(retinalIlluminanceImage,2), 3, 'uint8');
        theLMSexcitationsSequence = zeros(framesNum, size(retinalIlluminanceImage,1), size(retinalIlluminanceImage,2), 3, 'single');
        theSlices = zeros(framesNum, size(retinalIlluminanceImage,2), 'single');
       

        % Extract sequences
        for k = 1:framesNum
            theOI = theSequence.frameAtIndex(k);
            % Extract spatial support
            if (k == 1)
                spatialSupportMM = oiGet(theOI, 'spatial support', 'mm');
                optics = oiGet(theOI, 'optics');
                focalLength = opticsGet(optics, 'focal length');
                mmPerDegree = focalLength*tand(1)*1e3;
                spatialSupportDegsMatrix = spatialSupportMM/mmPerDegree;
                spatialSupportDegs.x = squeeze(spatialSupportDegsMatrix(1,:,1));
                spatialSupportDegs.y = squeeze(spatialSupportDegsMatrix(:,1,2));
                spatialSupportMicrons.x = 1e3*squeeze(spatialSupportMM(1,:,1));
                spatialSupportMicrons.y = 1e3*squeeze(spatialSupportMM(:,1,2));
            end
            retinalIrradianceImage = oiGet(theOI, 'energy');
            
            % Radiance to SRGB (scaled to max contrast)
            [~,~,retinalSRGBimageMaxContrast] = displayRadianceToDisplayRGB(retinalIrradianceImage, displaySPDs);
            theSRGBsequence(k,:,:,:) = uint8(retinalSRGBimageMaxContrast*255.0);
            
            % Compute retinal LMS excitations from the irradiance image
            rowsNum = size(retinalIrradianceImage,1);
            colsNum = size(retinalIrradianceImage,2);
            wavelenthsNum = size(retinalIrradianceImage,3);
            retinalPhotonImage = oiGet(theOI, 'photons');
            tmp = reshape(retinalPhotonImage, [rowsNum*colsNum wavelenthsNum]);
            tmp = tmp * coneQuantalEfficiencies;
            theLMSexcitationsSequence(k,:,:,:) = single(reshape(tmp, [rowsNum colsNum 3]));
            
            % retinal illuminance slice at mid-row
            retinalIlluminanceImage = oiGet(theOI, 'illuminance');
            theSlices(k,:) = single(squeeze(retinalIlluminanceImage(midRow,:)));
        end
        
    elseif (iscell(theSequence))
        % Dealing with a scene
        framesNum = numel(theSequence);
        luminanceImage = sceneGet(theSequence{1}, 'luminance');
        midRow = round(size(luminanceImage,1)/2);
        
        % Preallocate memory
        theSRGBsequence = zeros(framesNum , size(luminanceImage,1), size(luminanceImage,2), 3, 'uint8');
        theLMSexcitationsSequence = zeros(framesNum, size(luminanceImage,1), size(luminanceImage,2), 3, 'single');
        theSlices = zeros(framesNum, size(luminanceImage,2), 'single');
        
         % Extract sequences
        for k = 1:framesNum 
            theScene = theSequence{k};
            if (k == 1)
                % Extract spatial support
                spatialSupportMM = sceneGet(theScene, 'spatial support', 'mm');
                fovDegrees = sceneGet(theScene, 'wAngular');
                xSupport = squeeze(spatialSupportMM (1,:,1));
                spatialSupportDegs.x = xSupport / max(abs(xSupport(:))) * fovDegrees(1)/2;
                ySupport = squeeze(spatialSupportMM (:,1,2));
                spatialSupportDegs.y  = ySupport / max(abs(ySupport(:))) * fovDegrees(2)/2;
                spatialSupportMicrons.x = 1e3*squeeze(spatialSupportMM(1,:,1));
                spatialSupportMicrons.y = 1e3*squeeze(spatialSupportMM(:,1,2));
            end
            
            [~, sceneSRGBimage, ~, sceneLMSexcitationsImage] = ...
                sceneRepresentations(theScene, presentationDisplay);
            
            % SRGB sequence and LMSexcitations already available
            theSRGBsequence(k,:,:,:) = uint8(sceneSRGBimage*255.0);
            theLMSexcitationsSequence(k,:,:,:) = single(sceneLMSexcitationsImage);
            
            % luminance slice at mid-row
            luminanceImage = sceneGet(theScene, 'luminance');
            theSlices(k,:) = single(squeeze(luminanceImage(midRow,:)));
        end
         
    end
    
end
