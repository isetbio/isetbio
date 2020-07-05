function computeConeResponses(runParams, theConeMosaic, theOptics, ...
    recomputeNullResponses, ...
    instancesNum, stimDurationSeconds, stimulusFOVdegs, stimulusPixelsNum, ...
    spatialFrequenciesCPD, temporalFrequency, LMScontrast, ...
    saveDir)

    % Generate sine-wave scene as realized on a display
    stimColor = struct(...
        'backgroundChroma', [0.3, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', LMScontrast);
    
    stimTemporalParams = struct(...
        'temporalFrequencyHz', temporalFrequency, ...
        'stimDurationSeconds', stimDurationSeconds);
    
    % Determine phase increment at each stimulus frame so as to get required temporal frequency with the
    % constraint that the frame duration is evenly divisible by cone mosaic integration time
    [deltaPhaseDegs, stimTemporalParams.temporalFrequencyHz] = ...
        computeDeltaPhaseDegs(stimTemporalParams.temporalFrequencyHz, theConeMosaic.integrationTime);
    
    stimSpatialParams = struct(...
        'fovDegs', stimulusFOVdegs,...
        'pixelsNum', stimulusPixelsNum, ...
        'gaborPosDegs', [0 0], ...
        'gaborSpatialFrequencyCPD', 0, ...
        'gaborSigmaDegs', Inf, ... %stimulusFOVdegs/(2*4), ...%Inf, ...
        'gaborOrientationDegs', 0, ...
        'deltaPhaseDegs', deltaPhaseDegs);
    
    if (recomputeNullResponses)
        % Compute the null (zero contrast, no-noise responses)
        nullColor = stimColor;
        nullColor.lmsContrast = [0.0 0.0 0.0];
        
        % Generate the null stimulus frame scene sequence
        wavelengthSampling = theConeMosaic.pigment.wave;
        [theNullSceneFrames, presentationDisplay] = generateStimulusFrames(nullColor, stimSpatialParams, wavelengthSampling);
        
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
        % Extract the RGB sequence of optical images for the  null stimulus
        theStimulusRGBsequence = extractRGBsequence(theOIsequence, presentationDisplay);
        
        % Extract the stimulus and response time axes
        stimulusTimeAxis = theOIsequence.timeAxis;
        responseTimeAxis = theConeMosaic.timeAxis;
        
        % Save the null responses, and the photocurrent impulse responses/mean currents
        save(fullfile(saveDir, nullResponseFilename(runParams)), ...
             'theStimulusRGBsequence', ...
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
            LMScontrast(1), LMScontrast(2), LMScontrast(3), ...
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
        % Extract the RGB sequence of optical images for the  null stimulus
        [theStimulusRGBsequence, retinalIlluminanceSlice] = extractRGBsequence(theOIsequence, presentationDisplay);
       
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
        
        save(fullfile(saveDir,sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast), stimSpatialParams.gaborSpatialFrequencyCPD)), ...
            'theStimulusRGBsequence', ...
            'isomerizations', ...
            'photocurrents', ...
            '-v7.3');
        
        fprintf('Done in %2.1f minutes!\n', toc/60);
        clear('theStimulusRGBsequence');
        clear('isomerizations');
        clear('photocurrents');
        clear('theOIsequence');
        
    end % sfIndex

end

function [theRGBsequence,retinalIlluminanceSlice] = extractRGBsequence(theOIsequence, presentationDisplay)
    
    displaySPDs = displayGet(presentationDisplay, 'spd'); 
    theOI = theOIsequence.frameAtIndex(1);
    retinalIlluminanceImage = oiGet(theOI, 'illuminance');
    framesNum = theOIsequence.length;
    theRGBsequence = zeros(framesNum , size(retinalIlluminanceImage,1), size(retinalIlluminanceImage,2), 3, 'uint8');
    retinalIlluminanceSlice = zeros(framesNum, size(retinalIlluminanceImage,2), 'single');
    midRow = round(size(retinalIlluminanceImage,1)/2);
   
    for k = 1:framesNum
        theOI = theOIsequence.frameAtIndex(k);
        retinalIrradianceImage = oiGet(theOI, 'energy');
        retinalIlluminanceImage = oiGet(theOI, 'illuminance');
        retinalIlluminanceSlice(k,:) = single(squeeze(retinalIlluminanceImage(midRow,:)));
        [~, retinalSRGBimage] = displayRadianceToDisplayRGB(retinalIrradianceImage, displaySPDs);
        theRGBsequence(k,:,:,:) = uint8(retinalSRGBimage*255.0);
    end
    
end
