%
% RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction
%
%
function spatialTransferFunction(theMRGCMosaic, theOI, ...
    STFparamsStruct, theInputConeMosaicSTFResponsesFullFileName, ...
    theMRGCMosaicSTFResponsesFullFileName, varargin)

  p = inputParser;
  p.addParameter('computeInputConeMosaicResponses', false, @islogical);
  p.addParameter('computeMRGCMosaicResponses', false, @islogical);
  p.addParameter('debugInputConeMosaicPcurrentResponse', false, @islogical);
  p.addParameter('visualizeResponse', false, @islogical);
  p.addParameter('visualizeStimulusSequence', false, @islogical);
  p.addParameter('mRGCNonLinearityParams', [], @(x)(isempty(x))||(isstruct(x)));
  p.addParameter('customTemporalFrequencyAndContrast', [], @(x)(isempty(x))||(isstruct(x)));
  p.addParameter('validateScenes', false, @islogical);
  p.addParameter('visualizeCustomConeFundamentals', false, @islogical);
  

  % Execute the parser
  p.parse(varargin{:});
  computeInputConeMosaicResponses = p.Results.computeInputConeMosaicResponses;
  computeMRGCMosaicResponses = p.Results.computeMRGCMosaicResponses;
  mRGCNonLinearityParams = p.Results.mRGCNonLinearityParams;
  customTemporalFrequencyAndContrast = p.Results.customTemporalFrequencyAndContrast;
  visualizeResponse = p.Results.visualizeResponse;
  visualizeStimulusSequence = p.Results.visualizeStimulusSequence;
  debugInputConeMosaicPcurrentResponse = p.Results.debugInputConeMosaicPcurrentResponse;
  validateScenes = p.Results.validateScenes;
  visualizeCustomConeFundamentals = p.Results.visualizeCustomConeFundamentals;

  if (computeInputConeMosaicResponses)
    computePhotocurrent = false;
    pCurrentTemporalResolutionSeconds = [];
    osBiophysicalModelWarmUpTimeSeconds = [];
    osBiophysicalModelTimeStep = [];

    if (~isempty(mRGCNonLinearityParams)) && (strcmp(mRGCNonLinearityParams.type, 'photocurrent'))
          computePhotocurrent = true;
          osBiophysicalModelWarmUpTimeSeconds = mRGCNonLinearityParams.osBiophysicalModelWarmUpTimeSeconds;
          osBiophysicalModelTimeStep = mRGCNonLinearityParams.osBiophysicalModelTemporalResolutionSeconds;
          pCurrentTemporalResolutionSeconds = mRGCNonLinearityParams.pCurrentTemporalResolutionSeconds;
    end

    computeInputConeMosaicSTF(theMRGCMosaic, theOI, STFparamsStruct, ...
         theInputConeMosaicSTFResponsesFullFileName, ...
         computePhotocurrent, pCurrentTemporalResolutionSeconds, osBiophysicalModelWarmUpTimeSeconds, osBiophysicalModelTimeStep, ...
         customTemporalFrequencyAndContrast, ...
         visualizeResponse, visualizeStimulusSequence, debugInputConeMosaicPcurrentResponse, ...
         validateScenes, visualizeCustomConeFundamentals);
  end

  if (computeMRGCMosaicResponses)
    computeMRGCMosaicSTF(theMRGCMosaic, theInputConeMosaicSTFResponsesFullFileName, ...
      theMRGCMosaicSTFResponsesFullFileName, visualizeResponse, ...
      mRGCNonLinearityParams);
  end

end



function computeInputConeMosaicSTF(theMRGCMosaic, theOI, STFparamsStruct, theInputConeMosaicSTFResponsesFullFileName, ...
  computePhotocurrent,  pCurrentTemporalResolutionSeconds, photocurrentModelWarmUpTimeSeconds, osTimeStep, ...
  customTemporalFrequencyAndContrast, ...
  visualizeResponse, visualizeStimulusSequence, debugInputConeMosaicPcurrentResponse, ...
  validateScenes, visualizeCustomConeFundamentals)

  fprintf('Input cone mosaic STF responses will be saved in: \n%s\n', theInputConeMosaicSTFResponsesFullFileName);
  % Determine spatial frequencies examined
  minPixelsPerSpatialPeriod = 6;

  % Add the maxSF given the stimulus resolution, but only if this is <= 120CPD
  maxSF = min([120 1/(minPixelsPerSpatialPeriod*STFparamsStruct.resolutionDegs)]);

  spatialFrequenciesExamined = STFparamsStruct.sfSupport(STFparamsStruct.sfSupport <= maxSF);
  if (max(spatialFrequenciesExamined) < maxSF - 0.01)
    spatialFrequenciesExamined(numel(spatialFrequenciesExamined)+1) = maxSF;
  end

  % Orientations examined
  orientationsExamined = 0:STFparamsStruct.orientationDeltaDegs:(180-STFparamsStruct.orientationDeltaDegs);
  
  % Compute cone contrasts for desired chromaticity
  coneContrasts = ...
      visualStimulusGenerator.coneContrastsFromChromaticity(STFparamsStruct.chromaticity);


  % The defaults TF and achromatic contrast
  if (isempty(customTemporalFrequencyAndContrast))
    temporalFrequencyHz = 1.0;
    totalContrast = 0.75;
    backgroundLuminanceMultiplier = 1.0;
  else
    temporalFrequencyHz = customTemporalFrequencyAndContrast.temporalFrequencyHz;
    totalContrast = customTemporalFrequencyAndContrast.totalContrast;
    backgroundLuminanceMultiplier = customTemporalFrequencyAndContrast.backgroundLuminanceMultiplier;
  end



  stimParams = struct(...
    'backgroundChromaticity', STFparamsStruct.backgroundChromaticity, ...
    'backgroundLuminanceCdM2', STFparamsStruct.backgroundLuminanceCdM2, ...
    'contrast', totalContrast, ...
    'coneContrasts', coneContrasts, ...
    'sizeDegs', STFparamsStruct.sizeDegs, ...
    'positionDegs', STFparamsStruct.positionDegs, ...
    'resolutionDegs', STFparamsStruct.resolutionDegs, ...
    'spatialPhaseIncrementDegs', STFparamsStruct.spatialPhaseIncrementDegs, ...
    'temporalFrequencyHz', temporalFrequencyHz, ...
    'durationSeconds', 1.0/temporalFrequencyHz, ...
    'orientationDegs', orientationsExamined(1), ...
    'spatialFrequencyCPD', spatialFrequenciesExamined(1), ...
    'coneMosaicModulationBasedResponse',  true ...
    );

  % Generate presentation display, 20% luminance headroom
  viewingDistanceMeters = 4;

  thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
            theMRGCMosaic.inputConeMosaic.wave, ...
            STFparamsStruct.resolutionDegs, ...
            viewingDistanceMeters, ...
            'bitDepth', 20, ...
            'meanLuminanceCdPerM2', STFparamsStruct.backgroundLuminanceCdM2, ...
            'luminanceHeadroom', STFparamsStruct.displayLuminanceHeadroomPercentage);


  % Allocate memory (Single precision responses) to store all the cone mosaic responses
  [~, ~, spatialPhasesDegs, stimulusFrameSequenceTemporalSupportSeconds] = ...
      visualStimulusGenerator.driftingGratingModulationPatterns(stimParams);

  theInputConeMosaicSTFresponses = zeros(...
        numel(orientationsExamined), ...
        numel(spatialFrequenciesExamined), ...
        numel(spatialPhasesDegs), ...
        theMRGCMosaic.inputConeMosaic.conesNum, ...
        'single');

  if (STFparamsStruct.coneFundamentalsOptimizedForStimPosition)
      % Compute custom cone fundamentals
      maxConesNumForAveraging = 3;
      customConeFundamentals = visualStimulusGenerator.coneFundamentalsForPositionWithinConeMosaic(...
            theMRGCMosaic.inputConeMosaic, theOI, stimParams.positionDegs, stimParams.sizeDegs, maxConesNumForAveraging);
  end

  for iOri = 1:numel(orientationsExamined)
  for iSF = 1:numel(spatialFrequenciesExamined)
    
    % Get stim params
    stimParams.orientationDegs = orientationsExamined(iOri);
    stimParams.spatialFrequencyCPD = spatialFrequenciesExamined(iSF);

    % Feedback
    fprintf('Computing input cone mosaic STF for ORI = %d degrees and SF = %2.3f c/deg\n', ...
      stimParams.orientationDegs, stimParams.spatialFrequencyCPD);

    % Generate the spatial modulation patterns for all spatial phases of the drifting grating
    [theDriftingGratingSpatialModulationPatterns, spatialSupportDegs, spatialPhasesDegs, ...
      temporalSupportSeconds, temporalRamp] = visualStimulusGenerator.driftingGratingModulationPatterns(stimParams);

    if (STFparamsStruct.coneFundamentalsOptimizedForStimPosition)
        if ((iOri==1)&&(iSF==1))
            % Generate scenes for the different frames of the drifting grating and for the null stimulus
            [theDriftingGratingFrameScenes, theNullStimulusScene, ~, theConeFundamentalsStruct] = visualStimulusGenerator.stimulusFramesScenes(...
                  thePresentationDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
                  'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                  'customConeFundamentals', customConeFundamentals, ...
                  'announceEmployedConeFundamentals', true, ...
                  'validateScenes', validateScenes, ...
                  'visualizeCustomConeFundamentals', visualizeCustomConeFundamentals);
        else
            [theDriftingGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
                  thePresentationDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
                  'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                  'withPreviouslyComputedConeFundamentalsStruct', theConeFundamentalsStruct, ...
                  'announceEmployedConeFundamentals', true, ...
                  'validateScenes', ~true);
        end

    else
         % Generate scenes for the different frames of the drifting grating and for the null stimulus
        [theDriftingGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
                  thePresentationDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
                  'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                  'validateScenes', ~true);
    end

    % Compute input cone mosaic response to this orientation & spatial frequency
    stimulusPosition = 'mosaic-centered';

    % Set the mosaic integration time equal to the duration of one stimulus frame
    theMRGCMosaic.inputConeMosaic.integrationTime = temporalSupportSeconds(2)-temporalSupportSeconds(1);


    % Compute cone mosaic responses to each stimulus frame
    [theInputConeMosaicSTFresponses(iOri, iSF,:,:), theConeMosaicNullResponse] = ...
      RGCMosaicConstructor.helper.simulateExperiment.inputConeMosaicResponseToStimulusFrameSequence(...
        theMRGCMosaic, theOI, theNullStimulusScene, theDriftingGratingFrameScenes, ...
        stimulusPosition, stimParams.coneMosaicModulationBasedResponse, ...
        'visualizeResponse', visualizeResponse, ...
        'visualizeStimulusSequence', visualizeStimulusSequence, ...
        'thePresentationDisplayForVisualizingOpticalSceneOrImage', thePresentationDisplay, ...
        'stimulusInfoString', sprintf('ORI:%d degs, SF:%2.2f c/deg', stimParams.orientationDegs, stimParams.spatialFrequencyCPD));


    % Apply backgroundLuminanceMultiplier to the background cone excitations
    assert(backgroundLuminanceMultiplier == 1, sprintf('Why is the background luminance multiplier %f, not 1', backgroundLuminanceMultiplier));
    %theConeMosaicNullResponse = theConeMosaicNullResponse * backgroundLuminanceMultiplier;

    % Plot the LMS cone modulation responses
    if (debugInputConeMosaicPcurrentResponse)
      theLconeModulations = squeeze(theInputConeMosaicSTFresponses(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.lConeIndices));
      theMconeModulations = squeeze(theInputConeMosaicSTFresponses(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.mConeIndices));
      theSconeModulations = squeeze(theInputConeMosaicSTFresponses(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.sConeIndices));

      theLconeBackgroundExcitations = squeeze(theConeMosaicNullResponse(1,1,theMRGCMosaic.inputConeMosaic.lConeIndices));
      theMconeBackgroundExcitations = squeeze(theConeMosaicNullResponse(1,1,theMRGCMosaic.inputConeMosaic.mConeIndices));
      theSconeBackgroundExcitations = squeeze(theConeMosaicNullResponse(1,1,theMRGCMosaic.inputConeMosaic.sConeIndices));

      maxConeModulations = [max(abs(theLconeModulations(:))) max(abs(theMconeModulations(:))) max(abs(theSconeModulations(:))) ];
      meanConeExcitations = [mean(theLconeBackgroundExcitations(:)) mean(theMconeBackgroundExcitations(:)) mean(theSconeBackgroundExcitations(:))];

      figure(33);
      subplot(3,1,1)
      plot(temporalSupportSeconds, theLconeModulations, 'r-')
      title('Lcone response modulations');
      subplot(3,1,2)
      plot(temporalSupportSeconds, theMconeModulations, 'g-')
      title('Mcone response modulations');
      subplot(3,1,3)
      plot(temporalSupportSeconds, theSconeModulations, 'b-')
      title('Scone response modulations');

    end


    if (computePhotocurrent)
      stimulusPeriodDuration = 1/temporalFrequencyHz;
      nWarmUpPeriods = ceil(photocurrentModelWarmUpTimeSeconds/stimulusPeriodDuration);

      theConeMosaicResponseSequence = squeeze(theInputConeMosaicSTFresponses(iOri, iSF,:,:));

      % Transform cone modulations to cone excitations
      if (stimParams.coneMosaicModulationBasedResponse)
        theConeMosaicNullResponse = squeeze(theConeMosaicNullResponse);
        theConeMosaicNullResponse = reshape(theConeMosaicNullResponse, [1 numel(theConeMosaicNullResponse)]);

        % Transform from modulations to excitations
        %Rmod = (Recx-Ro)/Ro -> (Rmod+1)*Ro = RExc
        theConeMosaicResponseSequence = theConeMosaicResponseSequence  + 1;
        theConeMosaicResponseSequence = bsxfun(@times, theConeMosaicResponseSequence, theConeMosaicNullResponse);
      end

      % Compute the photocurrent response
      [theInputConeMosaicPhotocurrentTemporalSupportSeconds, ...
       theInputConeMosaicPhotocurrents(iOri, iSF,:,:), ...
       theInputConeMosaicBackgroundPhotocurrents(iOri, iSF,:)] = computePhotocurrents(...
            theMRGCMosaic.eccentricityDegs, temporalSupportSeconds, ...
            theConeMosaicResponseSequence, nWarmUpPeriods, ...
            pCurrentTemporalResolutionSeconds, osTimeStep, ...
            theMRGCMosaic.inputConeMosaic.coneTypes, ...
            debugInputConeMosaicPcurrentResponse);


      if (debugInputConeMosaicPcurrentResponse)
        theLconePhotocurrents = squeeze(theInputConeMosaicPhotocurrents(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.lConeIndices));
        theMconePhotocurrents = squeeze(theInputConeMosaicPhotocurrents(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.mConeIndices));
        theSconePhotocurrents = squeeze(theInputConeMosaicPhotocurrents(iOri, iSF,:,theMRGCMosaic.inputConeMosaic.sConeIndices));

        figure(34);
        subplot(3,1,1)
        plot(theInputConeMosaicPhotocurrentTemporalSupportSeconds, theLconePhotocurrents)
        subplot(3,1,2)
        plot(theInputConeMosaicPhotocurrentTemporalSupportSeconds, theMconePhotocurrents)
        subplot(3,1,2)
        plot(theInputConeMosaicPhotocurrentTemporalSupportSeconds, theSconePhotocurrents);
        drawnow;

      end

    end


   if (debugInputConeMosaicPcurrentResponse)
    pause
   end


  end % iSF
  end % iORI

  % Save results
  fprintf('Saving computed input cone mosaic STF responses to %s.\n', theInputConeMosaicSTFResponsesFullFileName);
  stimParams.orientationDegs = orientationsExamined;
  stimParams.spatialFrequencyCPD = spatialFrequenciesExamined;
  stimParams.spatialPhasesDegs = spatialPhasesDegs;
  stimParams.temporalSupportSeconds = temporalSupportSeconds;
  stimParams.temporalRamp = temporalRamp;

  save(theInputConeMosaicSTFResponsesFullFileName, ...
        'theMRGCMosaic', ...
        'STFparamsStruct', ...
        'stimParams', ...
        'theInputConeMosaicSTFresponses', ...
        '-v7.3');

  if (computePhotocurrent)
    % Append photocurrents
    theInputConeMosaicPhotocurrents = single(theInputConeMosaicPhotocurrents);
    theInputConeMosaicBackgroundPhotocurrents = single(theInputConeMosaicBackgroundPhotocurrents);

    save(theInputConeMosaicSTFResponsesFullFileName, ...
        'theInputConeMosaicPhotocurrentTemporalSupportSeconds', ...
        'theInputConeMosaicPhotocurrents', ...
        'theInputConeMosaicBackgroundPhotocurrents', ...
        '-append');
  end

end


function [temporalSupportPhotocurrent, theConePhotocurrents, theConeBackgroundPhotocurrents] = ...
          computePhotocurrents(eccentricityDegs, temporalSupportSeconds, theConeMosaicExcitationResponseSequence, ...
            nWarmUpPeriods, pCurrentTemporalResolutionSeconds, osTimeStep, coneTypes, debugInputConeMosaicPcurrentResponse)


  temporalSupportSecondsPeriodic = temporalSupportSeconds(1:end-1);
  temporalSupportSecondsPeriodic = reshape(temporalSupportSecondsPeriodic, [1 numel(temporalSupportSecondsPeriodic)]);

  nTimeBins = numel(temporalSupportSecondsPeriodic);
  for k = 1:nWarmUpPeriods
    temporalSupportSecondsPeriodic = cat(2, temporalSupportSecondsPeriodic, temporalSupportSeconds(1:nTimeBins)+k*temporalSupportSeconds(nTimeBins+1));
  end

  % Cone mosaic integration time
  coneMosaicIntegrationTime = temporalSupportSeconds(2)-temporalSupportSeconds(1);

  % Preallocate memory
  theSingleConeExcitations = theConeMosaicExcitationResponseSequence(1:end-1,1);
  temporalSupportPhotocurrent = computeSingleConePhotocurrentResponse( ...
      eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, nWarmUpPeriods, ...
      coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, false);

  theConePhotocurrents = zeros(numel(temporalSupportPhotocurrent), size(theConeMosaicExcitationResponseSequence,2));
  theConeBackgroundPhotocurrents = zeros(1, size(theConeMosaicExcitationResponseSequence,2));

  if (debugInputConeMosaicPcurrentResponse)
    for iConeIndex = 1:size(theConeMosaicExcitationResponseSequence,2)

      % Make the cone excitation response periodic by concatenating nWarmUpPeriods
      theSingleConeExcitations = theConeMosaicExcitationResponseSequence(1:end-1,iConeIndex);

      [~, theConePhotoCurrentDifferentialResponse, ...
       theConeBackgroundPhotoCurrent, ...
       theConeExcitationsSingleConePeriodic, ...
       photocurrentResponseTimeAxisPeriodic, ...
       thePcurrentResponsePeriodic, ...
       thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
                eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
                nWarmUpPeriods, coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, true);

      theConePhotocurrents(:, iConeIndex) = theConePhotoCurrentDifferentialResponse;
      theConeBackgroundPhotocurrents(iConeIndex) = theConeBackgroundPhotoCurrent;

      switch (coneTypes(iConeIndex))
        case cMosaic.LCONE_ID
          theColor = [231 187 192]/255;
        case cMosaic.MCONE_ID
          theColor = [128 207 218]/255;
        case cMosaic.SCONE_ID
          theColor = [0 0 1];
      end

      theResponse = squeeze(theConeMosaicExcitationResponseSequence(:,iConeIndex));
      m1 = max(theResponse(:));
      m2 = min(theResponse(:));
      coneContrast = (m1-m2)/(m1+m2);
      dt = (temporalSupportSeconds(2)-temporalSupportSeconds(1))/2;

      coneExcitationsRange = [0 1800];
      coneExcitationsTicks = 0:100:2000;

      coneExcitationsRange = [0 500];
      coneExcitationsTicks = 0:50:2000;

      coneExcitationRateRange = [0 20000];
      coneExcitationRateTicks = 0:2000:30000;

      coneExcitationRateTickLabels = cell(1, numel(coneExcitationRateTicks));
      for i = 1:numel(coneExcitationRateTicks)
        coneExcitationRateTickLabels{i} = sprintf('%2.0fk', coneExcitationRateTicks(i)/1000);
      end

      photocurrentRange = [-75 -15];
      deltaPhotocurrentRange = [-25 25];

      %photocurrentRange = [-80 -0];
      %deltaPhotocurrentRange = [-45 35]

      photocurrentTicks = -80:5:80;
      timeTicks = 0:100:10000;

      dTSeconds = temporalSupportSeconds(2)-temporalSupportSeconds(1);

      hFig = figure(33); clf;
      set(hFig, 'Position', [10 10 1600 1150], 'Color', [1 1 1]);
      subplot(2,3,3);

      stairs(temporalSupportSeconds*1e3, theResponse/dTSeconds, 'Color', 'k', 'LineWidth', 3.0);
      hold on;
      stairs(temporalSupportSeconds*1e3, theResponse/dTSeconds, 'Color', theColor, 'LineWidth', 1.5);
      plot([temporalSupportSeconds(1) temporalSupportSeconds(end)]*1e3, 0.5*(m1+m2)/dTSeconds*[1 1], '--', 'LineWidth', 1.0, 'Color', 'k', 'LineWidth', 1.0);
      ylabel('cone excitation rate (R*/sec)')
      xlabel('time (msec)');
      xtickangle(90);

      title(sprintf('response modulation: %2.2f', coneContrast), 'FontWeight', 'Normal');
      set(gca, 'YLim', coneExcitationRateRange, 'YTick', coneExcitationRateTicks, 'YTickLabel', coneExcitationRateTickLabels, ...
        'XTick', timeTicks, 'XLim', [temporalSupportSeconds(1) temporalSupportSeconds(end)+dt]*1e3, 'FontSize', 20);
      grid on; box off

      subplot(2,3,[1 2]);
      stairs(temporalSupportSecondsPeriodic*1e3, theConeExcitationsSingleConePeriodic/dTSeconds, 'Color', 'k', 'LineWidth', 3);
      hold on;
      stairs(temporalSupportSecondsPeriodic*1e3, theConeExcitationsSingleConePeriodic/dTSeconds, 'Color', theColor, 'LineWidth', 1.5);
      plot([temporalSupportSecondsPeriodic(1) temporalSupportSecondsPeriodic(end)]*1e3, 0.5*(m1+m2)/dTSeconds*[1 1], '--', 'LineWidth', 1.0, 'Color', 'k', 'LineWidth', 1.0);
      ylabel('cone excitation rate (R*/sec)')
      xlabel('time (msec)');
      title(sprintf('response modulation: %2.2f', coneContrast), 'FontWeight', 'Normal');
      xtickangle(90);
      set(gca, 'YLim', coneExcitationRateRange, 'YTick', coneExcitationRateTicks, 'YTickLabel', coneExcitationRateTickLabels, ...
        'XTick', timeTicks, 'XLim', [temporalSupportSecondsPeriodic(1) temporalSupportSecondsPeriodic(end)]*1e3, 'FontSize', 20);
      grid on; box off

      subplot(2,3,6);
      plot(temporalSupportPhotocurrent*1e3, theConePhotoCurrentDifferentialResponse, '-', 'Color', 'k', 'LineWidth', 3);
      hold on
      plot(temporalSupportPhotocurrent*1e3, theConePhotoCurrentDifferentialResponse, '-', 'Color', theColor, 'LineWidth', 1.5);
      plot(temporalSupportPhotocurrent*1e3,  temporalSupportPhotocurrent*0, 'k--');
      title(sprintf('(+): %2.2f, (-): %2.2f', max(theConePhotoCurrentDifferentialResponse), min(theConePhotoCurrentDifferentialResponse)), 'FontWeight', 'Normal');
      ylabel('differential pCurrent (pAmps)')
      xlabel('time (msec)');
      xtickangle(90);
      grid on; box off
      set(gca, 'YLim', deltaPhotocurrentRange, 'YTick', photocurrentTicks, 'XTick', timeTicks, 'XLim', [temporalSupportPhotocurrent(1) temporalSupportPhotocurrent(end)+dt]*1e3, 'FontSize', 20);


      subplot(2,3,[4 5]);
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentResponsePeriodic+theConeBackgroundPhotoCurrent, '-', 'Color', 'k', 'LineWidth', 3);
      hold on;
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentResponsePeriodic+theConeBackgroundPhotoCurrent, '-', 'Color', theColor, 'LineWidth', 1.5);
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentBackgroundResponseTransient, '--', 'Color', 'k', 'LineWidth', 1.0);
      set(gca, 'YLim', photocurrentRange , 'YTick', photocurrentTicks,  'XTick', timeTicks, 'XLim', [photocurrentResponseTimeAxisPeriodic(1) photocurrentResponseTimeAxisPeriodic(end)]*1e3, 'FontSize', 20);
      grid on; box off
      xtickangle(90);
      title(sprintf('(+): %2.2f, (-): %2.2f', max(theConePhotoCurrentDifferentialResponse), min(theConePhotoCurrentDifferentialResponse)), 'FontWeight', 'Normal');
      ylabel('pCurrent (pAmps)')
      xlabel('time (msec)');

      disp('Hit enter to continue')
      pause
    end % for iConeIndex
  else
    parfor iConeIndex = 1:size(theConeMosaicExcitationResponseSequence,2)
      % Make the cone excitation response periodic by concatenating nWarmUpPeriods
      theSingleConeExcitations = theConeMosaicExcitationResponseSequence(1:end-1,iConeIndex);

      [~, theConePhotoCurrentDifferentialResponse, theConeBackgroundPhotoCurrent] = ...
          computeSingleConePhotocurrentResponse(eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, nWarmUpPeriods, ...
            coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, true);

      theConePhotocurrents(:, iConeIndex) = theConePhotoCurrentDifferentialResponse;
      theConeBackgroundPhotocurrents(iConeIndex) = theConeBackgroundPhotoCurrent;
    end % parfor iConeIndex
  end

end


function [temporalSupportPhotocurrent, theConePhotoCurrentDifferentialResponse, backgroundPhotocurrent, ...
    theConeExcitationsPeriodic, photocurrentResponseTimeAxis, ...
    thePcurrentDifferentialPeriodicResponse, thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
            eccentricityDegs, theConeExcitations, temporalSupportSeconds, ...
            nWarmUpPeriods, coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, skipAssertions)

    % Make the cone excitation response periodic by concatenating nWarmUpPeriods
    theConeExcitationsPeriodic = reshape(theConeExcitations, [1 numel(theConeExcitations)]);
    for k = 1:nWarmUpPeriods
      theConeExcitationsPeriodic = cat(2, theConeExcitationsPeriodic, theConeExcitationsPeriodic(1:numel(theConeExcitations)));
    end

    % Convert excitation counts to excitation rates
    theConeExcitationsRatePeriodic = theConeExcitationsPeriodic(:) / coneMosaicIntegrationTime;
    backgroundConeExcitationRate = mean(theConeExcitationsRatePeriodic);

    % Compute the pCurrent response to the period stimulus
    [thePcurrentDifferentialPeriodicResponse, photocurrentResponseTimeAxis, thePcurrentBackgroundResponseTransient] = ...
        cMosaic.photocurrentFromConeExcitationRateUsingBiophysicalOSmodel(...
            sqrt(sum(eccentricityDegs(:).^2)), ...
            theConeExcitationsRatePeriodic, ...
            backgroundConeExcitationRate, ...
            coneMosaicIntegrationTime, ...
            pCurrentTemporalResolutionSeconds, ...
            'osTimeStepSeconds', osTimeStep, ...
            'skipAssertions', skipAssertions);

    backgroundPhotocurrent = thePcurrentBackgroundResponseTransient(end);

    % Only keep pCurrent response during the last stimulus period
    dToriginal = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    dT = photocurrentResponseTimeAxis(2)-photocurrentResponseTimeAxis(1);
    tOneStimulusCycle = temporalSupportSeconds(end)-temporalSupportSeconds(1);
    idx = find(photocurrentResponseTimeAxis >= photocurrentResponseTimeAxis(end)-(tOneStimulusCycle+0.5*dToriginal-dT));

    theConePhotoCurrentDifferentialResponse = thePcurrentDifferentialPeriodicResponse(idx);
    temporalSupportPhotocurrent = photocurrentResponseTimeAxis(idx);

    % Time support starts at 0 msec
    temporalSupportPhotocurrent = temporalSupportPhotocurrent - temporalSupportPhotocurrent(1);
end




function computeMRGCMosaicSTF(thePassedMRGCMosaic, theInputConeMosaicSTFResponsesFullFileName, ...
  theMRGCMosaicSTFResponsesFullFileName, visualizeComputedMRCMosaicSpatioTemporalResponse, mRGCNonLinearityParams)


  if (~isempty(mRGCNonLinearityParams)) && (strcmp(mRGCNonLinearityParams.type, 'photocurrent'))
    % Load photocurrent responses
    load(theInputConeMosaicSTFResponsesFullFileName, 'theMRGCMosaic', 'stimParams', ...
      'theInputConeMosaicSTFresponses', ...
      'theInputConeMosaicPhotocurrentTemporalSupportSeconds', ...
      'theInputConeMosaicPhotocurrents');

     % Check that all photocurrent response dimensions (other than time) agree with the cone excitation responses
     assert(size(theInputConeMosaicSTFresponses,1) == size(theInputConeMosaicPhotocurrents,1), ...
       'Mismatch in # of orientations');

     assert(size(theInputConeMosaicSTFresponses,2) == size(theInputConeMosaicPhotocurrents,2), ...
        'Mismatch in # of spatial frequencies');

     assert(size(theInputConeMosaicSTFresponses,4) == size(theInputConeMosaicPhotocurrents,4), ...
        'Mismatch in # of cones');

     theInputConeMosaicSTFresponses = theInputConeMosaicPhotocurrents;
     theConeMosaicResponseTemporalSupportSeconds = theInputConeMosaicPhotocurrentTemporalSupportSeconds;
  else

    fprintf('Loading computed input cone mosaic STF responses from %s.\n', theInputConeMosaicSTFResponsesFullFileName);
    load(theInputConeMosaicSTFResponsesFullFileName, 'theMRGCMosaic', 'stimParams', 'theInputConeMosaicSTFresponses');

    theConeMosaicResponseTemporalSupportSeconds = stimParams.temporalSupportSeconds;
  end


  % Assert that the input cone mosaic in the stored mRGC mosaic is the same as the input cone mosaic of the passed mRGC Mosaic
  assert(thePassedMRGCMosaic.inputConeMosaic.conesNum == theMRGCMosaic.inputConeMosaic.conesNum, ...
    'Number of cones in the loaded input cone mosaic STF responses file is not identicial to the that of input cone mosaic of the passed mRGC mosaic');
  assert(all(thePassedMRGCMosaic.inputConeMosaic.coneTypes == theMRGCMosaic.inputConeMosaic.coneTypes), ...
    'Cone types in the loaded input cone mosaic STF responses file are not identicial to those of the input cone mosaic of the passed mRGC mosaic');

  % Input cone mosaics are identical so switch to the passed mRGC mosaic which may have different characteristics of the mRGCs but identical
  % input cone mosaic
  clear 'theMRGCMosaic'
  theMRGCMosaic = thePassedMRGCMosaic;


  debugConnectivity = false;
  if (debugConnectivity)
      hFig = figure(44);
      ax = subplot(1,1,1);
      for theRGCindex = 1:theMRGCMosaic.rgcsNum

        centerConeIndices = theMRGCMosaic.singleCellConnectivityStats( theRGCindex, 'center', ...
            'inputConeIndicesOnly', true, ...
            'minConeWeightIncluded', 0.0001);

        surroundConeIndices = theMRGCMosaic.singleCellConnectivityStats( theRGCindex, 'surround', ...
            'inputConeIndicesOnly', true, ...
            'minConeWeightIncluded', 0.0001);

        plot(ax,...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,1), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,2), ...
            'bo');
        
        hold(ax, 'on')
        plot(ax,...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,1), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,2), ...
            'r.');
        drawnow;
        hold(ax, 'off')
      end % for theRGCindex 
  end % if (debugConnectivity)




  if (visualizeComputedMRCMosaicSpatioTemporalResponse)
    hFig = figure(55);
    clf;
    ax = subplot(1,1,1);
  else
      hFig = [];
      ax = [];
  end
  

  % Allocate memory
  coneMosaicResponseSize = [1 numel(theConeMosaicResponseTemporalSupportSeconds) theMRGCMosaic.inputConeMosaic.conesNum];
  theNoiseFreeConeMosaicSpatioTemporalExcitationsResponse = ...
            reshape(squeeze(theInputConeMosaicSTFresponses(1, 1, :,:)), coneMosaicResponseSize);

  [~, ~, theMRGCMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
                theNoiseFreeConeMosaicSpatioTemporalExcitationsResponse, theConeMosaicResponseTemporalSupportSeconds, ...
                'nonLinearityParams', mRGCNonLinearityParams);

  theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF = zeros(...
          numel(stimParams.orientationDegs), numel(stimParams.spatialFrequencyCPD), ...
          numel(theMRGCMosaicResponseTemporalSupportSeconds), theMRGCMosaic.rgcsNum);


  for iOri = 1:numel(stimParams.orientationDegs)  
    parfor iSF = 1:numel(stimParams.spatialFrequencyCPD)
        
        theNoiseFreeConeMosaicSpatioTemporalExcitationsResponse = ...
            reshape(squeeze(theInputConeMosaicSTFresponses(iOri, iSF, :,:)), coneMosaicResponseSize);

        [theNoiseFreeSpatioTemporalMRCMosaicResponse, ~, ...
            theMRGCMosaicResponseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
                theNoiseFreeConeMosaicSpatioTemporalExcitationsResponse, theConeMosaicResponseTemporalSupportSeconds, ...
                'nonLinearityParams', mRGCNonLinearityParams);

        % Normalize for contrast
        theNoiseFreeSpatioTemporalMRCMosaicResponse = theNoiseFreeSpatioTemporalMRCMosaicResponse /stimParams.contrast;
        mRGCMosaicActivationRange = [-1 1];

        if (visualizeComputedMRCMosaicSpatioTemporalResponse)
            for iTimeBin = 1:numel(theMRGCMosaicResponseTemporalSupportSeconds)
              theFrameResponse = theNoiseFreeSpatioTemporalMRCMosaicResponse(1,iTimeBin,:);
              theMRGCMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'activation', theFrameResponse, ...
                'activationRange', mRGCMosaicActivationRange);
              drawnow;
            end
        end

      theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF(iOri, iSF, :,:) = theNoiseFreeSpatioTemporalMRCMosaicResponse;
    end % for iSF
  end % for iOri

  save(theMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic', 'stimParams', 'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF', 'theMRGCMosaicResponseTemporalSupportSeconds');
end
