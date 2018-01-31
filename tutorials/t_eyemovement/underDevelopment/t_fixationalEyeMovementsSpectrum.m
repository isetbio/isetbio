function t_fixationalEyeMovementsSpectrum

    emDurationSeconds = 1; sampleTimeSeconds = 1/1000; nTrials = 100;
    
    % Initialize object
    fixEMobj = fixationalEM();
    
    microSaccadeType = 'none'; % 'heatmap/fixation based';
    
    configNames = {...
%        'zeroFeedback' ...
%         'doubleFeedback'
%        'zeroFeedbackOneThirdGamma' ...
%        'positionalNoiseOnly' ...
        'defaultParams' ...
%         'lowerPersistence' ...
%         'higherPersistence' ...
%         'higherPersistence_less_feedback' ...
%         'higherPersistence_more_feedback' ...
%         'lowerPersistence_less_feedback' ...
%         'lowerPersistence_more_feedback' ...   
        };
    
    velocityMeasurementIntervalSeconds = 41/1000;
    for configNo = 1:numel(configNames)
        % Set params for positional noise only
        fixEMobj = paramsForConfig(configNames{configNo}, microSaccadeType, fixEMobj);

        % Compute the random seed to get reproducible results
        %fixEMobj.randomSeed = 3457;
        % Compute the emPaths
        fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials);

        % Visualize
        figNo = configNo;
        hFig = analyzeAndDisplayResults(figNo, fixEMobj, configNames{configNo}, velocityMeasurementIntervalSeconds);
        
        % Export to PDF
        NicePlot.exportFigToPNG(sprintf('%s.png', configNames{configNo}), hFig, 300);
    end
end

function fixEMobj = paramsForConfig(configName, microSaccadeType, fixEMobj)
    fixEMobj.setDefaultParams();
    fixEMobj.stabilizationSeconds = 3;
    fixEMobj.microSaccadeType = microSaccadeType;
    
    switch configName
        case 'defaultParams' 
   
        case 'higherPersistence'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * 0.8;
            
        case 'lowerPersistence'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * (1/0.8);
            
        case 'higherPersistence_less_feedback'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * 0.8;
            fixEMobj.feedbackGain = fixEMobj.feedbackGain * 0.8;
            
        case 'higherPersistence_more_feedback'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * 0.8;
            fixEMobj.feedbackGain = fixEMobj.feedbackGain * (1/0.8);
           
        case 'lowerPersistence_less_feedback'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * (1/0.8);
            fixEMobj.feedbackGain = fixEMobj.feedbackGain * 0.8;
            
        case 'lowerPersistence_more_feedback'
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlGamma = fixEMobj.controlGamma * (1/0.8);
            fixEMobj.feedbackGain = fixEMobj.feedbackGain * (1/0.8);
            
        case 'positionalNoiseOnly' 
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.controlNoiseStd = 0;
            
        case 'zeroFeedback' 
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.feedbackGain = 0;
        
        case 'doubleFeedback' 
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.feedbackGain = fixEMobj.feedbackGain * (1.5);
            
        case 'zeroFeedbackOneThirdGamma' 
            fixEMobj.microSaccadeType = 'none';
            fixEMobj.feedbackGain = 0;
            fixEMobj.controlGamma = fixEMobj.controlGamma / 3;
            
        otherwise
            error('Unknown config name: ''%s''.', configName);
    end
end

function hFig = analyzeAndDisplayResults(figNo, fixEMobj, configName,  velocityMeasurementIntervalSeconds)

    % Analyze the emPaths
    nTrials = size(fixEMobj.emPathsArcMin,1);
    timeStepSeconds = fixEMobj.timeAxis(2)-fixEMobj.timeAxis(1);
    
    emPosRange = [-35 35];
    emPosDelta = 1;
    
    for trialNo = 1:nTrials
        emPathArcMin = squeeze(fixEMobj.emPathsArcMin(trialNo,:,:));
        xPosArcMin = squeeze(emPathArcMin(:,1));
        yPosArcMin = squeeze(emPathArcMin(:,2));
        
        % Perform displacement analysis
        [dispX, ~, scrambledIntervalsDispX, DrandomWalkX, timeLagsMilliseconds] = ...
            fixEMobj.performDisplacementAnalysis(xPosArcMin', fixEMobj.timeAxis, 'mode', 'D1');
        [dispY, ~, scrambledIntervalsDispY, DrandomWalkY, ~] = ...
            fixEMobj.performDisplacementAnalysis(yPosArcMin', fixEMobj.timeAxis, 'mode', 'D1');
        
        [dispX2, ~, scrambledIntervalsDispX2, DrandomWalkX2, ~] = ...    
            fixEMobj.performDisplacementAnalysis(xPosArcMin', fixEMobj.timeAxis, 'mode', 'D2');
        [dispY2, ~, scrambledIntervalsDispY2, DrandomWalkY2, ~] = ...    
            fixEMobj.performDisplacementAnalysis(yPosArcMin', fixEMobj.timeAxis, 'mode', 'D2');
        
        % Compute velocity
        velocityArcMinPerSecond = fixEMobj.computeVelocity(emPathArcMin', velocityMeasurementIntervalSeconds, timeStepSeconds);

        % Perform spectral analysis
        samplingRate = 1/timeStepSeconds;
        [pSDX,frequencyAxis] = periodogram(xPosArcMin,[],[], samplingRate);
        [pSDY,~] = periodogram(yPosArcMin,[],[], samplingRate);
    
        if (trialNo == 1)
            powerSpectralDensityX = zeros(nTrials, length(pSDX));
            powerSpectralDensityY = zeros(nTrials, length(pSDX));
            displacementX = zeros(nTrials, length(dispX));
            displacementY = zeros(nTrials, length(dispY));
            scrambledIntervalsDisplacementX = zeros(nTrials, length(scrambledIntervalsDispX));
            scrambledIntervalsDisplacementY = zeros(nTrials, length(scrambledIntervalsDispY));
            displacementX2 = zeros(nTrials, length(dispX2));
            displacementY2 = zeros(nTrials, length(dispY2));
            scrambledIntervalsDisplacementX2 = zeros(nTrials, length(scrambledIntervalsDispX2));
            scrambledIntervalsDisplacementY2 = zeros(nTrials, length(scrambledIntervalsDispY2));
        end 
         
        % Accumulate over trials
        displacementX(trialNo,:) = dispX;
        displacementY(trialNo,:) = dispY;
        scrambledIntervalsDisplacementX(trialNo,:) = scrambledIntervalsDispX;
        scrambledIntervalsDisplacementY(trialNo,:) = scrambledIntervalsDispY;
        
        displacementX2(trialNo,:) = dispX2;
        displacementY2(trialNo,:) = dispY2;
        scrambledIntervalsDisplacementX2(trialNo,:) = scrambledIntervalsDispX2;
        scrambledIntervalsDisplacementY2(trialNo,:) = scrambledIntervalsDispY2;
        
        powerSpectralDensityX(trialNo,:) = pSDX;
        powerSpectralDensityY(trialNo,:) = pSDY;
    end % trialNo
    
    % Mean over trials
    powerSpectralDensityX = squeeze(mean(powerSpectralDensityX,1));
    powerSpectralDensityY = squeeze(mean(powerSpectralDensityY,1));
    displacementX = squeeze(mean(displacementX,1));
    displacementY = squeeze(mean(displacementY,1));
    %displacement = sqrt(displacementX.^2 + displacementY.^2);
    scrambledIntervalsDisplacementX = squeeze(mean(scrambledIntervalsDisplacementX,1));
    %scrambledIntervalsDisplacementY = squeeze(mean(scrambledIntervalsDisplacementY,1));

    displacementX2 = squeeze(mean(displacementX2,1));
    displacementY2 = squeeze(mean(displacementY2,1));
    scrambledIntervalsDisplacementX2 = squeeze(mean(scrambledIntervalsDisplacementX2,1));
    scrambledIntervalsDisplacementY2 = squeeze(mean(scrambledIntervalsDisplacementY2,1));
    displacementX2degs = displacementX2/(60^2);
    displacementY2degs = displacementY2/(60^2);
    
    scrambledIntervalsDisplacementX2degs = scrambledIntervalsDisplacementX2/(60^2);
    scrambledIntervalsDisplacementY2degs = scrambledIntervalsDisplacementY2/(60^2);
    displacement2Degs = displacementX2degs + displacementY2degs;
    scrambledIntervalsDisplacement2Degs = scrambledIntervalsDisplacementX2degs + scrambledIntervalsDisplacementY2degs;
     
    % Plot the results
    hFig = figure(figNo); clf; 
    set(hFig, 'Position', [0 0 1150 950], 'Name', configName, 'Color', [1 1 1]);

    % emPlots - position (last trial)
    subplot(3,3,[1 2]);
    plot(fixEMobj.timeAxis, xPosArcMin, 'r-', 'LineWidth', 1.5); hold on;
    plot(fixEMobj.timeAxis, yPosArcMin, 'b-', 'LineWidth', 1.5); 
    plot(fixEMobj.timeAxis, yPosArcMin*0, 'k-');
    hold off;
    grid on; 
    timeTicks = 0:0.2:100;
    set(gca, 'XLim', [fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        'XTick', timeTicks, 'XTickLabel', {}, ...
        'YLim', emPosRange, 'FontSize', 14);
    ylabel('position (arc min)');
    legend({'x-pos', 'y-pos'})
    title('x/y eye position (last trial)');
    
    % Velocity (last trial)
    subplot(3,3,[4 5]);
    plot(fixEMobj.timeAxis, velocityArcMinPerSecond, 'b-', 'LineWidth', 1.5);
    hold on;
    subjectDR.mean = 30;
    subjectDR.std = 17;
    plot([fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], (subjectDR.mean-subjectDR.std)*[1 1], 'k--', 'LineWidth', 1.0);
    subjectDG.mean = 89;
    subjectDG.std = 63;
    plot([fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], (subjectDG.mean+subjectDG.std)*[1 1], 'k--', 'LineWidth', 1.0);
    hold off;
    grid on
    legend({'model', 'Cherici et al: subjectDR', 'Cherici et al: subjectDG'}, 'Location', 'NorthEast')
    set(gca, 'XLim', [fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        'XTick', timeTicks, 'XTickLabel', sprintf('%2.2f\n', timeTicks), ...
        'YLim', [0 60*10], 'YTick', [0:100:60*100], 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('velocity (arc min / sec)');
    title(sprintf('eye velocity (last trial, velocity interval: %2.0f msec)', velocityMeasurementIntervalSeconds*1000));
    
    subplot(3,3,3);
    % Last trial
    plot(xPosArcMin, yPosArcMin, 'k-', 'LineWidth', 1.0); hold on
    plot([0 0], emPosRange, 'k-'); plot(emPosRange, [0 0], 'k-');  hold off;
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, 'XTick', [-100:10:100], 'YTick', [-100:10:100], 'XTickLabel', {}, 'FontSize', 14);
    grid on; box on
    axis 'square'
    ylabel('y-position (arc min)');
    title('eye movement path (last trial)');
    
    subplot(3,3,6);
    hold on
    % fixation span
    cmap = jet(1024);
    [emLikelihoodMap, emLikelihoodMapSupportX, emLikelihoodMapSupportY] = ...
        fixEMobj.computeLikelihoodMap(fixEMobj.emPathsArcMin, emPosRange, emPosDelta);
    contourf(emLikelihoodMapSupportX, emLikelihoodMapSupportY, emLikelihoodMap, 0:0.05:1, 'LineColor', 'none'); hold on;
    plot([0 0], emPosRange, 'k-'); plot(emPosRange, [0 0], 'k-');  hold off;
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, 'XTick', [-100:10:100], 'YTick', [-100:10:100], 'FontSize', 14);
    grid on; box on
    axis 'square'
    colormap(cmap);
    hh = colorbar('north', 'Ticks', [0:0.2:1.0], 'TickLabels', {'.0', '.2', '.4', '.6', '.8', '1.'}, 'Color', [1 1 1]);
    xlabel('x-position (arc min)');
    ylabel('y-position (arc min)');
    title(sprintf('fixation span (%d trials)', nTrials));
    
    % Spectral density plots
    subplot(3,3,7);
    plot(frequencyAxis,10*log10(powerSpectralDensityX), 'r-', 'LineWidth', 1.5); hold on
    plot(frequencyAxis,10*log10(powerSpectralDensityY), 'b-', 'LineWidth', 1.5); hold off
    set(gca, 'XScale', 'log', 'XLim', [1 300], 'YLim', [-50 10], 'XTick', [1 3 10 30 100 300], 'YTick', [-100:10:100], 'FontSize', 14);
    grid on;
    axis 'square'
    xlabel('frequency (Hz)')
    ylabel('power density (dB)')
    title(sprintf('power spectral density\n(%d trials)', nTrials));
    
    % Displacement^2 plot
    subplot(3,3,8);
    plot(timeLagsMilliseconds, displacement2Degs, 'k-', 'LineWidth', 1.5); hold on;
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacement2Degs, 'k--', 'LineWidth', 1.5); hold off;
    grid on
    set(gca, 'XLim', [1 1000], 'YLim', [1e-4 1e-1], 'XTick', [1 3 10 30 100 300 1000], ...
        'YTick', [1e-4 1e-3 1e-2 1e-1], 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    axis 'square'
    xlabel('time interval (ms)')
    ylabel('mean displacement^2 (deg^2)');
    legend({'emPath', 'scrambled'}, 'Location', 'SouthEast')
    title(sprintf('mean squared displacement\n(%d trials)', nTrials));
    
    % Linear displacement plot
    subplot(3,3,9);
    plot(timeLagsMilliseconds, displacementX, 'r-', 'LineWidth', 2); hold on;
    plot(timeLagsMilliseconds, displacementY, 'b-', 'LineWidth', 2);
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacementX, 'k--', 'LineWidth', 1.5); hold off;
    grid on
    set(gca, 'XLim', [2 1000], 'YLim', [0.1 10], 'XTick', [3 10 30 100 300 1000], ...
        'YTick', [0.03 0.1 0.3 1 3 10], 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    axis 'square'
    xlabel('time interval (ms)')
    ylabel('mean displacement (min arc)')
    legend({'emPath (x)', 'emPath (y)', 'scrambled (x)'}, 'Location', 'SouthEast');
    title(sprintf('mean displacement\n(%d trials)', nTrials));
    drawnow;
end

