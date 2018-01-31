function t_fixationalEyeMovementsExploreModel

    emDurationSeconds = 3; sampleTimeSeconds = 1/1000; nTrials = 100;
    
    % Initialize object
    fixEMobj = fixationalEM();
    fixEMobj.microSaccadeType = 'none';
    feedbackSteepness = fixEMobj.feedbackSteepness * 1.0;
    
    defaults.gamma = fixEMobj.controlGamma;
    defaults.feedback = fixEMobj.feedbackGain;
    
    %controlGammaValues = [0.25]
    feedbackGainValues = [0 0.10 0.125  0.15  0.175 0.20 0.225];
    
    controlGammaValues = [0.15 0.17 0.20   0.25   0.5];
    %feedbackGainValues = [0.15];
    
    [feedbakGain, controlGamma] = meshgrid(feedbackGainValues, controlGammaValues);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', size(controlGamma,1), ...
                   'colsNum', size(controlGamma,2), ...
                   'heightMargin',   0.04, ...
                   'widthMargin',    0.01, ...
                   'leftMargin',     0.03, ...
                   'rightMargin',    0.001, ...
                   'bottomMargin',   0.04, ...
                   'topMargin',      0.03);
               
    hFig1 = figure(1); clf; set(hFig1, 'Position', [10 10 1900 1100], 'Color', [1 1 1]); 
    hFig2 = figure(2); clf; set(hFig2, 'Position', [10 10 1900 1100], 'Color', [1 1 1]); 
    hFig3 = figure(3); clf; set(hFig3, 'Position', [10 10 1900 1100], 'Color', [1 1 1]); 
    hFig4 = figure(4); clf; set(hFig4, 'Position', [10 10 1900 1100], 'Color', [1 1 1]); 
    hFig5 = figure(5); clf; set(hFig5, 'Position', [10 10 1900 1100], 'Color', [1 1 1]); 
    
    emPosRange = [-20 20];
    emPosDelta = 1;
    
    for iVar = 1:numel(controlGamma)
        % Set params for positional noise only
        fixEMobj.controlGamma = controlGamma(iVar);
        fixEMobj.feedbackGain = feedbakGain(iVar);
        fixEMobj.feedbackSteepness = feedbackSteepness;
        
        % Compute the random seed to get reproducible results
        fixEMobj.randomSeed = 3457;
        % Compute the emPaths
        fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials);

        % Analyze the results
        velocityMeasurementIntervalSeconds = 41/1000;
        d = analyzeResults(fixEMobj, velocityMeasurementIntervalSeconds, emPosRange, emPosDelta);
 
        plotEMpath(hFig1, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, fixEMobj.timeAxis, d.emPathArcMin, d.emLikelihoodMapSupportX, d.emLikelihoodMapSupportY, d.emLikelihoodMap);
        plotVelocity(hFig2, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, fixEMobj.timeAxis, d.velocityArcMinPerSecond);
        plotPowerSpectralDensity(hFig3, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, d.frequencyAxis, d.powerSpectralDensityX, d.powerSpectralDensityY);
        plotDisplacementD1(hFig4, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, d.timeLagsMilliseconds, d.displacementX, d.displacementY, d.scrambledIntervalsDisplacementX);
        plotDisplacementD2(hFig5, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, d.timeLagsMilliseconds, d.displacement2Degs, d.scrambledIntervalsDisplacement2Degs);
    end
    
    % Export to PDF
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'emPath', feedbackSteepness), hFig1, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'velocity', feedbackSteepness), hFig2, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'psd', feedbackSteepness), hFig3, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'd1', feedbackSteepness), hFig4, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'd2', feedbackSteepness), hFig5, 300);
end

function d = analyzeResults(fixEMobj, velocityMeasurementIntervalSeconds, emPosRange, emPosDelta)

    % Analyze the emPaths
    nTrials = size(fixEMobj.emPathsArcMin,1);
    timeStepSeconds = fixEMobj.timeAxis(2)-fixEMobj.timeAxis(1);
    
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
        velocity = fixEMobj.computeVelocity(emPathArcMin', velocityMeasurementIntervalSeconds, timeStepSeconds);

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
            velocityArcMinPerSecond = zeros(nTrials, length(velocity));
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
        
        velocityArcMinPerSecond(trialNo,:) = velocity;
    end % trialNo
    
    % Mean over trials
    d.velocityArcMinPerSecond = squeeze(mean(velocityArcMinPerSecond,1));
    
    d.powerSpectralDensityX = 10*log10(squeeze(mean(powerSpectralDensityX,1)));
    d.powerSpectralDensityY = 10*log10(squeeze(mean(powerSpectralDensityY,1)));
    d.frequencyAxis = frequencyAxis;
    
    d.displacementX = squeeze(mean(displacementX,1));
    d.displacementY = squeeze(mean(displacementY,1));
    %displacement = sqrt(displacementX.^2 + displacementY.^2);
    d.scrambledIntervalsDisplacementX = squeeze(mean(scrambledIntervalsDisplacementX,1));
    %scrambledIntervalsDisplacementY = squeeze(mean(scrambledIntervalsDisplacementY,1));

    displacementX2 = squeeze(mean(displacementX2,1));
    displacementY2 = squeeze(mean(displacementY2,1));
    scrambledIntervalsDisplacementX2 = squeeze(mean(scrambledIntervalsDisplacementX2,1));
    scrambledIntervalsDisplacementY2 = squeeze(mean(scrambledIntervalsDisplacementY2,1));
    displacementX2degs = displacementX2/(60^2);
    displacementY2degs = displacementY2/(60^2);
    
    scrambledIntervalsDisplacementX2degs = scrambledIntervalsDisplacementX2/(60^2);
    scrambledIntervalsDisplacementY2degs = scrambledIntervalsDisplacementY2/(60^2);
    d.displacement2Degs = displacementX2degs + displacementY2degs;
    d.scrambledIntervalsDisplacement2Degs = scrambledIntervalsDisplacementX2degs + scrambledIntervalsDisplacementY2degs;
    d.timeLagsMilliseconds = timeLagsMilliseconds;
    d.emPathArcMin = emPathArcMin;
    [d.emLikelihoodMap, d.emLikelihoodMapSupportX, d.emLikelihoodMapSupportY] = ...
        fixEMobj.computeLikelihoodMap(fixEMobj.emPathsArcMin, emPosRange, emPosDelta);
end


function plotEMpath(hFig, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, timeAxis, emPathArcMin, emLikelihoodMapSupportX, emLikelihoodMapSupportY, emLikelihoodMap)
    [row,col] = ind2sub(size(controlGamma), iVar);
    nRows = size(controlGamma,1);
    figure(hFig);
    subplot('Position', subplotPosVectors(row,col).v);
    emPosRange = [-10 10];
    emPosTicks = [-100:5:100];
    idx = find(timeAxis <= 0.2);
    imagesc(emLikelihoodMapSupportX, emLikelihoodMapSupportY, emLikelihoodMap, [0 1]); hold on;
    plot(emPathArcMin(idx,1), emPathArcMin(idx,2), 'r-', 'LineWidth', 1.0); 
    plot([0 0], emPosRange, 'k-'); plot(emPosRange, [0 0], 'k-');  hold off;
    colormap(gray(1024));
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, 'XTick', emPosTicks, 'YTick', emPosTicks, 'FontSize', 14);
    grid on; box on
    axis 'square'
    
    if (row == size(controlGamma,1))
        xlabel('x-position (arc min)');
    else
        set(gca, 'XTickLabel', {});
    end
    
    if (col == 1)
        ylabel('y-position (arc min)');
    else
        set(gca, 'YTickLabel', {});
    end
    
    if ((defaults.gamma == controlGamma(iVar)) && (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)), 'Color', [1 0 0]);
    else
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)));
    end
end

function plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)    
    
    [row,col] = ind2sub(size(controlGamma), iVar);
    nRows = size(controlGamma,1);
    figure(hFig);
    subplot('Position', subplotPosVectors(row,col).v);
    plot(timeAxis, velocityArcMinPerSecond, 'b-', 'LineWidth', 1.5);
    hold on;
    subjectDR.mean = 30;
    subjectDR.std = 17;
    plot([timeAxis(1) timeAxis(end)], (subjectDR.mean-subjectDR.std)*[1 1], 'k--', 'LineWidth', 1.0);
    subjectDG.mean = 89;
    subjectDG.std = 63;
    plot([timeAxis(1) timeAxis(end)], (subjectDG.mean+subjectDG.std)*[1 1], 'k--', 'LineWidth', 1.0);
    hold off;
    grid on
    timeTicks = 0:0.5:100;
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)], ...
        'XTick', timeTicks, 'XTickLabel', sprintf('%2.2f\n', timeTicks), ...
        'YLim', [0 60*10], 'YTick', [0:100:60*100], 'FontSize', 14);
    axis 'square'
    
    if (row == size(controlGamma,1))
        xlabel('time (sec)');
    else
        set(gca, 'XTickLabel', {});
    end
    
    if (col == 1)
        ylabel('vel. (arcmin/sec)');
    else
        set(gca, 'YTickLabel', {});
    end
    
    if ((defaults.gamma == controlGamma(iVar)) && (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)), 'Color', [1 0 0]);
    else
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)));
    end
    
end

function plotPowerSpectralDensity(hFig, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, frequencyAxis,powerSpectralDensityX,powerSpectralDensityY)   

    [row,col] = ind2sub(size(controlGamma), iVar);
    nRows = size(controlGamma,1);
    figure(hFig);
    subplot('Position', subplotPosVectors(row,col).v);
    plot(frequencyAxis,powerSpectralDensityX, 'r-', 'LineWidth', 1.5); hold on
    plot(frequencyAxis,powerSpectralDensityY, 'b-', 'LineWidth', 1.5); hold off
    set(gca, 'XScale', 'log', 'XLim', [1 300], 'YLim', [-50 10], 'XTick', [1 3 10 30 100 300], 'YTick', [-100:10:100], 'FontSize', 14);
    grid on;
    axis 'square'
    if (row == size(controlGamma,1))
        xlabel('freq. (Hz)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('power(dB)');
    else
        set(gca, 'YTickLabel', {});
    end
    
    if ((defaults.gamma == controlGamma(iVar)) && (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)), 'Color', [1 0 0]);
    else
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)));
    end
end

function plotDisplacementD2(hFig, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, timeLagsMilliseconds, displacement2Degs, scrambledIntervalsDisplacement2Degs)

    [row,col] = ind2sub(size(controlGamma), iVar);
    nRows = size(controlGamma,1);
    figure(hFig);
    subplot('Position', subplotPosVectors(row,col).v);
    
    plot(timeLagsMilliseconds, displacement2Degs, 'b-', 'LineWidth', 1.5); hold on;
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacement2Degs, 'k--', 'LineWidth', 1.5); hold off;
    grid on
    set(gca, 'XLim', [2 1000], 'YLim', [0.3*1e-4 1e-1], 'XTick', [3 10 30 100 300 1000], ...
        'YTick', [1e-4 1e-3 1e-2 1e-1], 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    axis 'square'
    if (row == size(controlGamma,1))
        xlabel('dt (ms)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('mean D^2 (deg^2)');
    else
        set(gca, 'YTickLabel', {});
    end
    
    if ((defaults.gamma == controlGamma(iVar)) && (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)), 'Color', [1 0 0]);
    else
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)));
    end
end

function plotDisplacementD1(hFig, subplotPosVectors, iVar, controlGamma, feedbakGain, defaults, timeLagsMilliseconds, displacementX, displacementY, scrambledIntervalsDisplacementX)

    [row,col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row,col).v);
    
    plot(timeLagsMilliseconds, displacementX, 'r-', 'LineWidth', 2); hold on;
    plot(timeLagsMilliseconds, displacementY, 'b-', 'LineWidth', 2);
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacementX, 'k--', 'LineWidth', 1.5); hold off;
    grid on
    set(gca, 'XLim', [0 500], 'YLim', [0 8], 'XTick', [0:100:1000], ...
        'YTick', [0:10], 'XScale', 'linear', 'YScale', 'linear', 'FontSize', 14);
    axis 'square'
    if (row == size(controlGamma,1))
        xlabel('dt (ms)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('mean D (min arc)');
    else
        set(gca, 'YTickLabel', {});
    end
    
    if ((defaults.gamma == controlGamma(iVar)) && (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)), 'Color', [1 0 0]);
    else
        title(sprintf('control:%0.2f, feedback:%0.3f', controlGamma(iVar), feedbakGain(iVar)));
    end
end

