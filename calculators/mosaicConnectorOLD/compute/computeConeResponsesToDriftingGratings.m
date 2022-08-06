function computeConeResponsesToDriftingGratings(runParams, ...
    stimColor,  stimTemporalParams, stimSpatialParams, ...
    theConeMosaic, theOptics, ...
    recomputeNullResponses, ...
    instancesNum, ...
    opticsPostFix, PolansSubjectID, ...
    saveDir, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('saveCornealStimulusSequence', true, @islogical);
    p.addParameter('saveRetinalStimulusSequence', true, @islogical);
    p.parse(varargin{:});
    saveCornealStimulusSequence = p.Results.saveCornealStimulusSequence;
    saveRetinalStimulusSequence = p.Results.saveRetinalStimulusSequence;
    
    % Determine phase increment at each stimulus frame so as to get required temporal frequency with the
    % constraint that the frame duration is evenly divisible by cone mosaic integration time
    [stimSpatialParams.deltaPhaseDegs, stimTemporalParams.temporalFrequencyHz] = ...
        computeDeltaPhaseDegs(stimTemporalParams.temporalFrequencyHz, theConeMosaic.integrationTime);
    
    % Whether to display plots of the generated scene and optical image
    debugStimulusGeneration = false;
    
    if (recomputeNullResponses)
        % Compute the null (zero contrast, no-noise responses)
        nullColor = stimColor;
        nullColor.lmsContrast = [0.0 0.0 0.0];
        
        % Generate the null stimulus frame scene sequence
        wavelengthSampling = theConeMosaic.pigment.wave;
        [theNullSceneFrames, presentationDisplay, stimSpatialParams.pixelsNum] = generateStimulusFrames(nullColor, stimSpatialParams, wavelengthSampling);
        
        if (saveCornealStimulusSequence)
            % Generate the stimulus SRGBsequence
            [theStimulusSRGBsequence, ~, theStimulusLMSexcitationsSequence, ...
                stimulusSpatialSupportDegs, stimulusSpatialSupportMicrons] = extractSRGBsequence(theNullSceneFrames, presentationDisplay, theConeMosaic.qe);
        end
        
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
        
        if (saveRetinalStimulusSequence)
            % Extract the RGB sequence of optical images and LMS excitation images
            [theOISRGBsequence, ~, theRetinalLMSexcitationsSequence, ...
                retinalSpatialSupportDegs, retinalSpatialSupportMicrons] = extractSRGBsequence(theOIsequence, presentationDisplay, theConeMosaic.qe);
        end
        
        % Extract the stimulus and response time axes
        stimulusTimeAxis = theOIsequence.timeAxis;
        responseTimeAxis = theConeMosaic.timeAxis;
        
        % Assemble the null response filename
        theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix, PolansSubjectID);
        
        % Save the null responses, and the photocurrent impulse responses/mean currents
        save(fullfile(saveDir, theNullResponseFileName), ...
             'stimColor', 'stimTemporalParams', 'stimSpatialParams', ...
             'saveCornealStimulusSequence', 'saveRetinalStimulusSequence', ...
             'isomerizationsNull', ...
             'photocurrentsNull', ...
             'responseTimeAxis', ...
             'stimulusTimeAxis', ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents', ...
             '-v7.3');
         
         % Append corneal stimulus sequence
         if (saveCornealStimulusSequence)
             save(fullfile(saveDir, theNullResponseFileName), ...
                'theStimulusSRGBsequence', ...
                'theStimulusLMSexcitationsSequence', ...
                'stimulusSpatialSupportDegs', 'stimulusSpatialSupportMicrons', ...
                '-append');
         end
         
         % Append retinal stimulus sequence
         if (saveRetinalStimulusSequence)
             save(fullfile(saveDir, theNullResponseFileName), ...
                'theOISRGBsequence', ...
                'theRetinalLMSexcitationsSequence', ...
                'retinalSpatialSupportDegs', 'retinalSpatialSupportMicrons', ...
                '-append');
         end
         
         fprintf('Done in %2.1f minutes\n', toc/60);
         
         % Clear some space from memory
         clear('theStimulusRGBsequence');
         clear('theStimulusLMSexcitationsSequence');
         clear('theOISRGBsequence');
         clear('theRetinalLMSexcitationsSequence');
         clear('isomerizationsNull');
         clear('photocurrentsNull');
         clear('theOIsequence')
    else
        fprintf('\nLoading null responses ...');
        tic
        % Load the photocurrent impulse responses and mean currents
        theNullResponseFileName = nullResponseFilename(runParams, opticsPostFix, PolansSubjectID);
        load(fullfile(saveDir, theNullResponseFileName), ...
             'osImpulseResponseFunctions', ...
             'osMeanCurrents');
         fprintf('Done in %2.1f minutes\n', toc/60);
    end
    
    

    % Extract tested spatial frequencies
    spatialFrequenciesCPD = stimSpatialParams.testedSpatialFrequenciesCPD;
    
    % Compute responses
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        % Test spatial frequency
        stimSpatialParams.gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Generate the test (spatial frequency) stimulus frame scene sequence
        wavelengthSampling = theConeMosaic.pigment.wave;
        [theSceneFrames, presentationDisplay] = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling);
        
        if (debugStimulusGeneration) || (saveCornealStimulusSequence)
            % Generate the stimulus SRGBsequence
            [theStimulusSRGBsequence, sceneLuminanceSlice, theStimulusLMSexcitationsSequence, ...
                stimulusSpatialSupportDegs, stimulusSpatialSupportMicrons] = extractSRGBsequence(theSceneFrames, presentationDisplay, theConeMosaic.qe);

            if (debugStimulusGeneration)
                sceneLuminanceRange = [0 100];
                maxSceneLuminance = updateStimulusDebugFigure(4000, sceneLuminanceSlice, sceneLuminanceRange, ...
                    spatialFrequenciesCPD, maxSceneLuminance, sfIndex, 'scene luminance');
            end
        end
        
        
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
        
        if (debugStimulusGeneration) || (saveRetinalStimulusSequence)
            % Extract the RGB sequence of optical images and LMS excitation images
            [theOISRGBsequence, retinalIlluminanceSlice, theRetinalLMSexcitationsSequence, ...
                retinalSpatialSupportDegs, retinalSpatialSupportMicrons] = extractSRGBsequence(theOIsequence, presentationDisplay, theConeMosaic.qe);

            if (debugStimulusGeneration)
                retinalIlluminanceRange = [1 2];
                maxRetinalIlluminance = updateStimulusDebugFigure(6000, retinalIlluminanceSlice, retinalIlluminanceRange, ...
                    spatialFrequenciesCPD, maxRetinalIlluminance, sfIndex, 'retinal illuminance');
            end
        end
        
        % Assemble the test response filename
        theTestResponseFileName = sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, stimColor.lmsContrast, opticsPostFix, PolansSubjectID), ...
            stimSpatialParams.gaborSpatialFrequencyCPD);
        
        save(fullfile(saveDir,theTestResponseFileName), ...
            'stimColor', 'stimTemporalParams', 'stimSpatialParams', ...
            'saveCornealStimulusSequence', 'saveRetinalStimulusSequence', ...
            'isomerizations', ...
            'photocurrents', ...
            '-v7.3');
        
        % Append corneal stimulus sequence
        if (saveCornealStimulusSequence)
            save(fullfile(saveDir, theTestResponseFileName), ...
                'theStimulusSRGBsequence', ...
                'theStimulusLMSexcitationsSequence', ...
                'stimulusSpatialSupportDegs', 'stimulusSpatialSupportMicrons', ...
                '-append');
        end
        
        % Append retinal stimulus sequence
        if (saveRetinalStimulusSequence)
            save(fullfile(saveDir, theTestResponseFileName), ...
                'theOISRGBsequence', ...
                'theRetinalLMSexcitationsSequence', ...
                'retinalSpatialSupportDegs', 'retinalSpatialSupportMicrons', ...
                '-append');
        end
         
        fprintf('Done in %2.1f minutes!\n', toc/60);
        clear('theStimulusSRGBsequence');
        clear('isomerizations');
        clear('photocurrents');
        clear('theOIsequence');
        
    end % sfIndex
end

function [theSRGBsequence, theSlices, theLMSexcitationsSequence, spatialSupportDegs, spatialSupportMicrons] = ...
    extractSRGBsequence(theSequence, presentationDisplay, coneQuantalEfficiencies)

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
                fovDegrees(1) = sceneGet(theScene, 'hfov');
                fovDegrees(2) = sceneGet(theScene, 'vfov');
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

function maxSlice = updateStimulusDebugFigure(figNo, theSlices, theSlicesRange, spatialFrequenciesCPD, maxSlice, sfIndex, yAxisLabel)
     figure(figNo);
     subplot(3,5,sfIndex)
     imagesc(theSlices);
     set(gca, 'CLim', theSlicesRange);
     maxSlice(sfIndex) = max(theSlices(:));

     subplot(3,5,15);
     plot(spatialFrequenciesCPD(1:sfIndex), maxSlice(1:sfIndex), 'ks-');
     ylabel(yAxisLabel);
     xlabel('c/deg');
     set(gca, 'XLim', [0.1 70], 'XScale', 'log');
     colormap(gray);
end

