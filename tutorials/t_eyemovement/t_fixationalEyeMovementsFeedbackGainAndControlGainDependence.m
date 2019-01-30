function t_fixationalEyeMovementsFeedbackGainAndControlGainDependence
% Explore how fixationalEM drift parameters affect emPath characteristics
%
% Syntax:
%   t_fixationalEyeMovementsFeedbackGainAndControlGainDependence
%
% Description:
%    Explore how the following key characteristics of emPaths
%       (i) velocity
%      (ii) fixation span
%     (iii) power spectal density
%      (iv) displacement analysis
%    depend on the following @fixationalEM drift component parameters
%     - the feedback gain parameter and
%     - the control gain parameter
%
%    The script generates eye movement paths with no micro-saccades to
%    probe the drift model in isolation. The 4 figures generated show how
%    one of the 4 key characteristics changes as the feedback gain/control
%    gain change. The panel with the red title corresponds to the default
%    values of feedback/control gain, (in Mergenthaler&Engbert, 2007).
%
%    Note that there are several combinations of feedback gain and control
%    gain, which result in drift velocities and fixation spans that are
%    well within the range reported by Cherici et al (2012) 
%    "Precision of sustained fixation in trained and untrained observers"
%
%    If one were interested to study how performance on a task depends on
%    eye movements one could use a number of combinations of feedback gain
%    and control gain to generate a population of subjects with slightly
%    different eye movement patterns. This tutorial provides a good start
%    for this.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History
%    02/06/18  npc  Wrote it.
%    02/07/18  npc  Comments.
%    07/12/18  jnm  Formatting. Written as a function, so formatted as one.

% Generate 500 5-second emPaths with a sample time of 1 msec
emDurationSeconds = 5;
sampleTimeSeconds = 1 / 1000;
nTrials = 500;

% Examined values of feedback and control gain
feedbackGainValues = [0, 0.10, 0.125, 0.15, 0.175, 0.20];
controlGammaValues = [0.17, 0.20, 0.25, 0.5];
[feedbakGainGrid, controlGammaGrid] = ...
    meshgrid(feedbackGainValues, controlGammaValues);

% Initialize object
fixEMobj = fixationalEM();
fixEMobj.microSaccadeType = 'none';

% Keep the default values of control and feedback
defaults.gamma = fixEMobj.controlGamma;
defaults.feedback = fixEMobj.feedbackGain;

% Subplot arrangement
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', size(controlGammaGrid, 1), ...
    'colsNum', size(controlGammaGrid, 2), ...
    'heightMargin', 0.02, ...
    'widthMargin', 0.01, ...
    'leftMargin', 0.035, ...
    'rightMargin', 0.001, ...
    'bottomMargin', 0.04, ...
    'topMargin', 0.02);

% Prepare figures
hFig1 = figure(1);
clf;
set(hFig1, 'Position', [10 10 1400 950], 'Color', [1 1 1]);
hFig2 = figure(2);
clf;
set(hFig2, 'Position', [100 10 1400 950], 'Color', [1 1 1]);
hFig3 = figure(3);
clf;
set(hFig3, 'Position', [200 10 1400 950], 'Color', [1 1 1]);
hFig4 = figure(4);
clf;
set(hFig4, 'Position', [300 10 1400 950], 'Color', [1 1 1]);
% hFig5 = figure(5);
% clf;
% set(hFig5, 'Position', [400 10 1400 950], 'Color', [1 1 1]);

for iVar = 1:numel(controlGammaGrid)
    % Set varied params
    fixEMobj.controlGamma = controlGammaGrid(iVar);
    fixEMobj.feedbackGain = feedbakGainGrid(iVar);
    fixEMobj.randomSeed = 3457;

    % Compute the emPaths
    computeVelocity = true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
        nTrials, computeVelocity, ...
        'useParfor', false);

    % Analyze the emPaths
    emPosRange = [-30 30];
    emPosDelta = 0.5;
    d = analyzeResults(fixEMobj, emPosRange, emPosDelta);

    % Plot different aspects of the emPaths
    % Plot the fixation span across all trials and a single emPath
    plotFixationSpan(hFig1, subplotPosVectors, iVar, controlGammaGrid, ...
        feedbakGainGrid, defaults, fixEMobj.timeAxis, emPosRange, ...
        d.emPathArcMin, d.fixationMap, d.fixationMapSupportX, ...
        d.fixationMapSupportY, d.fixationMapXSlice, d.fixationMapYSlice);

    % Plot the velocity
    plotVelocity(hFig2, subplotPosVectors, iVar, controlGammaGrid, ...
        feedbakGainGrid, defaults, fixEMobj.timeAxis, ...
        d.velocityArcMinPerSecond);

    % PLot the power spectral density
    plotPowerSpectralDensity(hFig3, subplotPosVectors, iVar, ...
        controlGammaGrid, feedbakGainGrid, defaults, d.frequencyAxis, ...
        d.powerSpectralDensityX, d.powerSpectralDensityY);

    % Plot the mean squared displacement (in deg^2)
    plotDisplacementD2(hFig4, subplotPosVectors, iVar, ...
        controlGammaGrid, feedbakGainGrid, defaults, ...
        d.timeLagsMilliseconds, d.displacement2Degs, ...
        d.scrambledIntervalsDisplacement2Degs);

    % Plot the mean displacement (in min arc)
    % plotDisplacementD1(hFig5, subplotPosVectors, iVar, controlGamma, ...
    %    feedbakGain, defaults, d.timeLagsMilliseconds, ...
    %    d.displacementX, d.displacementY, ...
    %    d.scrambledIntervalsDisplacementX);
end

% Export to PDF
exportFigsToPDF = false;
if (exportFigsToPDF)
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'emPath', ...
       feedbackSteepness), hFig1, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'velocity', ...
       feedbackSteepness), hFig2, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'psd', ...
       feedbackSteepness), hFig3, 300);
    NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'd1', ...
       feedbackSteepness), hFig4, 300);
%     NicePlot.exportFigToPDF(sprintf('%s_steepness_%2.2f.pdf', 'd2', ...
%        feedbackSteepness), hFig5, 300);
end

end

function d = analyzeResults(fixEMobj, emPosRange, emPosDelta)
% Perform displacement, fixationMap and spectral analysis and compute 
% the mean over all trials.
%
% Syntax:
%   d = analyzeResults(fixEMobj, emPosRange, emPosDelta)
%
% Description:
%    Analyze the emPaths.
%
% Inputs:
%    FixEMobj   - Object. A fixational eye movements object.
%    emPosRange - Vector. The eye movement position range.
%    emPosDelta - Numeric. The Delta for eye movement positions.
%
% Outputs:
%    d          - Struct. A structure containing the relevant results for
%                 the emPaths, including fixation map, and other attributes
%
% Optional key/value pairs:
%    None.
%
    nTrials = size(fixEMobj.emPosArcMin, 1);
    timeStepSeconds = fixEMobj.timeAxis(2) - fixEMobj.timeAxis(1);

    for trialNo = 1:nTrials
        emPathArcMin = squeeze(fixEMobj.emPosArcMin(trialNo, :, :));
        xPosArcMin = squeeze(emPathArcMin(:, 1));
        yPosArcMin = squeeze(emPathArcMin(:, 2));

        % Perform displacement analysis
        [dispX, ~, scrambledIntervalsDispX, DrandomWalkX, ...
            timeLagsMilliseconds] = fixEMobj.performDisplacementAnalysis(...
            xPosArcMin', fixEMobj.timeAxis, 'mode', 'D1');
        [dispY, ~, scrambledIntervalsDispY, DrandomWalkY, ~] = ...
            fixEMobj.performDisplacementAnalysis(yPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D1');

        [dispX2, ~, scrambledIntervalsDispX2, DrandomWalkX2, ~] = ...
            fixEMobj.performDisplacementAnalysis(xPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D2');
        [dispY2, ~, scrambledIntervalsDispY2, DrandomWalkY2, ~] = ...
            fixEMobj.performDisplacementAnalysis(yPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D2');

        % Perform spectral analysis
        samplingRate = 1 / timeStepSeconds;
        [pSDX, frequencyAxis] = periodogram(xPosArcMin, [], [], samplingRate);
        [pSDY, ~] = periodogram(yPosArcMin, [], [], samplingRate);

        if (trialNo == 1)
            powerSpectralDensityX = zeros(nTrials, length(pSDX));
            powerSpectralDensityY = zeros(nTrials, length(pSDX));
            displacementX = zeros(nTrials, length(dispX));
            displacementY = zeros(nTrials, length(dispY));
            scrambledIntervalsDisplacementX = ...
                zeros(nTrials, length(scrambledIntervalsDispX));
            scrambledIntervalsDisplacementY = ...
                zeros(nTrials, length(scrambledIntervalsDispY));
            displacementX2 = zeros(nTrials, length(dispX2));
            displacementY2 = zeros(nTrials, length(dispY2));
            scrambledIntervalsDisplacementX2 = ...
                zeros(nTrials, length(scrambledIntervalsDispX2));
            scrambledIntervalsDisplacementY2 = ...
                zeros(nTrials, length(scrambledIntervalsDispY2));
            % velocityArcMinPerSecond = zeros(nTrials, length(velocity));
        end

        % Accumulate over trials
        displacementX(trialNo, :) = dispX;
        displacementY(trialNo, :) = dispY;
        scrambledIntervalsDisplacementX(trialNo, :) = scrambledIntervalsDispX;
        scrambledIntervalsDisplacementY(trialNo, :) = scrambledIntervalsDispY;

        displacementX2(trialNo, :) = dispX2;
        displacementY2(trialNo, :) = dispY2;
        scrambledIntervalsDisplacementX2(trialNo, :) = ...
            scrambledIntervalsDispX2;
        scrambledIntervalsDisplacementY2(trialNo, :) = ...
            scrambledIntervalsDispY2;

        powerSpectralDensityX(trialNo, :) = pSDX;
        powerSpectralDensityY(trialNo, :) = pSDY;

        % velocityArcMinPerSecond(trialNo, :) = velocity;
    end % trialNo

    velocityArcMinPerSecond = fixEMobj.velocityArcMin;
    % fixEMobj.computeVelocity(emPathArcMin', ...
    %    velocityMeasurementIntervalSeconds, timeStepSeconds);

    % Mean over trials
    d.velocityArcMinPerSecond = squeeze(mean(velocityArcMinPerSecond, 1));

    d.powerSpectralDensityX = ...
        10 * log10(squeeze(mean(powerSpectralDensityX, 1)));
    d.powerSpectralDensityY = ...
        10 * log10(squeeze(mean(powerSpectralDensityY, 1)));
    d.frequencyAxis = frequencyAxis;

    d.displacementX = squeeze(mean(displacementX, 1));
    d.displacementY = squeeze(mean(displacementY, 1));
    % displacement = sqrt(displacementX .^ 2 + displacementY .^ 2);
    d.scrambledIntervalsDisplacementX = ...
        squeeze(mean(scrambledIntervalsDisplacementX, 1));
    % scrambledIntervalsDisplacementY = ...
    %    squeeze(mean(scrambledIntervalsDisplacementY, 1));

    displacementX2 = squeeze(mean(displacementX2, 1));
    displacementY2 = squeeze(mean(displacementY2, 1));
    scrambledIntervalsDisplacementX2 = ...
        squeeze(mean(scrambledIntervalsDisplacementX2, 1));
    scrambledIntervalsDisplacementY2 = ...
        squeeze(mean(scrambledIntervalsDisplacementY2, 1));
    displacementX2degs = displacementX2 / (60 ^ 2);
    displacementY2degs = displacementY2 / (60 ^ 2);

    scrambledIntervalsDisplacementX2degs = ...
        scrambledIntervalsDisplacementX2 / (60 ^ 2);
    scrambledIntervalsDisplacementY2degs = ...
        scrambledIntervalsDisplacementY2 / (60 ^ 2);
    d.displacement2Degs = displacementX2degs + displacementY2degs;
    d.scrambledIntervalsDisplacement2Degs = ...
        scrambledIntervalsDisplacementX2degs ...
        + scrambledIntervalsDisplacementY2degs;
    d.timeLagsMilliseconds = timeLagsMilliseconds;
    d.emPathArcMin = emPathArcMin;

    maxDurationSeconds = 1.0;
    [d.fixationMap, d.fixationMapSupportX, d.fixationMapSupportY, ...
        d.fixationMapXSlice, d.fixationMapYSlice] = ...
        fixEMobj.computeFixationMap(fixEMobj.timeAxis, ...
        fixEMobj.emPosArcMin, emPosRange, emPosDelta, ...
        'maxDurationSeconds', maxDurationSeconds);
end

function plotFixationSpan(hFig, subplotPosVectors, iVar, controlGamma, ...
    feedbakGain, defaults, timeAxis, emPosRange, emPathArcMin, ...
    fixationMap, fixationMapSupportX, fixationMapSupportY, ...
    fixationMapXSlice, fixationMapYSlice)
% Plot the fixation span across all trials and a single eye movement path
%
% Syntax:
%   plotFixationSpan(hFig, subplotPosVectors, iVar, controlGamma, ...
%       feedbackGain, defaults, timeAxis, emPosRange, emPathArcMin, ...
%       fixationMap, fixationMapSupportX, fixationMapSupportY, ...
%       fixationMapXSlice, fixationMapYSlice)
%
% Description:
%    Plot the eye movement path using the provided parameters.
%
% Inputs:
%    hFig                - Figure.  The figure handle.
%    subplotPosVectors   - Matrix. A matrix of structures, where each
%                          structure contains a vector.
%    iVar                - Numeric. The index variable for the plot.
%    controlGamma        - Matrix. The control gamma grid.
%    feedbackGain        - Matrix. The feedback gain grid.
%    defaults            - Struct. The defaults structure.
%    timeAxis            - Vector. A time axis vector.
%    emPosRange          - Vector. A 1x2 matrix containing the minimum and
%                          maximum values for the eye movement position.
%    emPathArcMin        - Matrix. A matrix containing the emPath in arcmin
%    fixationMap         - Matrix. A matrix containing the fixation map.
%    fixationMapSupportX - Vector. A vector containing the X axis support
%                          for the fixation map.
%    fixationMapSupportY - Vector. A vector containing the Y axis support
%                          for the fixation map.
%    fixationMapXSlice   - Vector. A vector containing a X slice of the
%                          fixation map.
%    fixationMapYSlice)  - Vector. A vector containing a Y slice of the
%                          fixation map.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    [row, col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row, col).v);
    emPosTicks = [-100:10:100];
    % Plot the fixation map
    contourf(fixationMapSupportX, fixationMapSupportY, fixationMap, ...
        0:0.05:1, 'LineColor', [.5 0.5 0.5]);
    hold on;
    plot([0 0], emPosRange, 'k-');
    plot(emPosRange, [0 0], 'k-');
    plot(fixationMapSupportX, emPosRange(1) + fixationMapXSlice * ...
        emPosRange(2) * 0.9, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
    plot(emPosRange(1) + fixationMapYSlice * emPosRange(2) * 0.9...
        , fixationMapSupportY, '-', 'Color', [0 0 1], 'LineWidth', 1.5);

    % Plot 1 second of a path
    idx = find(timeAxis <= 1.0);
    plot(emPathArcMin(idx, 1), emPathArcMin(idx, 2), 'g-', 'LineWidth', 1.5);
    hold off;
    colormap(brewermap(1024, 'Greys'));
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, 'XTick', emPosTicks, ...
        'YTick', emPosTicks, 'FontSize', 12);
    grid on;
    box on
    axis 'square';
    axis 'xy';

    if (row == size(controlGamma, 1))
        xlabel('x-position (arc min)');
    else
        set(gca, 'XTickLabel', {});
    end

    if (col == 1)
        ylabel('y-position (arc min)');
    else
        set(gca, 'YTickLabel', {});
    end

    if ((defaults.gamma == controlGamma(iVar)) && ...
            (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'Color', [1 0 0], 'FontSize', 10);
    else
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'FontSize', 10);
    end

    drawnow;
end

function plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, ...
    feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)
% Plot the mean eye movement velocity across all trials along with the
% velocity upper and lower bounds velocities from a pool of subjects from
% Cherici et al.
%
% Syntax:
%   plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, ...
%       feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)
%
% Description:
%    Plot the velocity figure.
%
% Inputs:
%    hFig                    - Figure. The figure handle.
%    subplotPosVectors       - Matrix. A matrix of structures, where each
%                              structure contains a vector.
%    iVar                    - Numeric. The index variable for the plot.
%    controlGamma            - Matrix. The control gamma grid.
%    feedbackGain            - Matrix. The feedback gain grid.
%    defaults                - Struct. The defaults structure.
%    timeAxis                - Vector. A time axis vector.
%    velocityArcMinPerSecond - Vector. a vector containing the velocities
%                              in arcmin per second.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    [row, col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row, col).v);
    plot(timeAxis, velocityArcMinPerSecond, 'k-', 'LineWidth', 1.5);
    hold on;
    
    % The mean and std for subject DR are from Table 3 of Cherici et al
    subjectDR.mean = 30;
    subjectDR.std = 17;
    
    plot([timeAxis(1) timeAxis(end)], ...
        (subjectDR.mean - subjectDR.std) * [1 1], 'k--', 'LineWidth', 1.0);
    
    % The mean and std for subject DG are from Table 3 of Cherici et al
    subjectDG.mean = 89;
    subjectDG.std = 63;
    
    plot([timeAxis(1) timeAxis(end)], ...
        (subjectDG.mean + subjectDG.std) * [1 1], 'k:', 'LineWidth', 1.0);
    hold off;
    grid on
    timeTicks = 0:0.5:100;
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)], ...
        'XTick', timeTicks, 'XTickLabel', sprintf('%2.2f\n', timeTicks), ...
        'YLim', [0 60 * 10], 'YTick', [0:100:60 * 100], 'FontSize', 12);
    axis 'square'

    if (row == size(controlGamma, 1))
        xlabel('time (sec)');
    else
        set(gca, 'XTickLabel', {});
    end

    if (col == 1)
        ylabel('velocity (arcmin/sec)');
    else
        set(gca, 'YTickLabel', {});
    end

    legend({'model', 'Cherici et al (best subj.)', ...
        'Cherici et al (worst subj.)'}, 'Location', 'NorthWest');
    if ((defaults.gamma == controlGamma(iVar)) && ...
            (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'Color', [1 0 0], 'FontSize', 10);
    else
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'FontSize', 10);
    end
    drawnow;
end

function plotPowerSpectralDensity(hFig, subplotPosVectors, iVar, ...
    controlGamma, feedbakGain, defaults, frequencyAxis, ...
    powerSpectralDensityX, powerSpectralDensityY)
% Plot the power spectral density for X and Y components of fixational eye
% movemetns.
%
% Syntax:
%   plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, ...
%       feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)
%
% Description:
%    Plot the velocity figure.
%
% Inputs:
%    hFig                  - Figure. The figure handle.
%    subplotPosVectors     - Matrix. A matrix of structures, where each
%                            structure contains a vector.
%    iVar                  - Numeric. The index variable for the plot.
%    controlGamma          - Matrix. The control gamma grid.
%    feedbackGain          - Matrix. The feedback gain grid.
%    defaults              - Struct. The defaults structure.
%    frequencyAxis         - Vector. A frequency axis vector.
%    powerSpectralDensityX - Vector. A vector containing the X axis support
%                            for the spectral density map.
%    powerSpectralDensityY - Vector. A vector containing the Y axis support
%                            for the spectral density map.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    [row, col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row, col).v);
    plot(frequencyAxis, powerSpectralDensityX, 'r-', 'LineWidth', 1.5);
    hold on
    plot(frequencyAxis, powerSpectralDensityY, 'b-', 'LineWidth', 1.5);
    hold off
    set(gca, 'XScale', 'log', 'XLim', [1 300], 'YLim', [-50 10], ...
        'XTick', [1 3 10 30 100 300], 'YTick', [-100:10:100], 'FontSize', 12);
    grid on;
    axis 'square'
    if (row == size(controlGamma, 1))
        xlabel('freq. (Hz)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('power(dB)');
    else
        set(gca, 'YTickLabel', {});
    end

    legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');

    if ((defaults.gamma == controlGamma(iVar)) && ...
            (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'Color', [1 0 0], 'FontSize', 10);
    else
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'FontSize', 10);
    end
    drawnow;
end

function plotDisplacementD2(hFig, subplotPosVectors, iVar, ...
    controlGamma, feedbakGain, defaults, timeLagsMilliseconds, ...
    displacement2Degs, scrambledIntervalsDisplacement2Degs)
% Plot the mean displacement squared as a function of time lag
%
% Syntax:
%   plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, ...
%       feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)
%
% Description:
%    Plot the velocity figure.
%
% Inputs:
%    hFig                                - Figure. The figure handle.
%    subplotPosVectors                   - Matrix. A matrix of structures,
%                                          wherein each structure contains
%                                          a vector for subplot positioning
%    iVar                                - Numeric. The index variable for
%                                          the plot.
%    controlGamma                        - Matrix. The control gamma grid.
%    feedbackGain                        - Matrix. The feedback gain grid.
%    defaults                            - Struct. The defaults structure.
%    timeLagsMilliseconds                - Vector. A vector containing the
%                                          lag in milliseconds.
%    displacement2Degs                   - Vector. A vector containing the
%                                          displacement converted into
%                                          units of degrees.
%    scrambledIntervalsDisplacement2Degs - Vector. A vector containing the
%                                          scrambled interval displacements
%                                          converted into units of degrees.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    [row, col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row, col).v);

    plot(timeLagsMilliseconds, displacement2Degs, 'k-', 'LineWidth', 1.5);
    hold on;
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacement2Degs, 'k--', ...
        'LineWidth', 1.5);
    hold off;
    grid on
    set(gca, 'XLim', [2 1000], 'YLim', [0.3 * 1e-4 1e-1], ...
        'XTick', [3 10 30 100 300 1000], 'YTick', [1e-4 1e-3 1e-2 1e-1], ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
    axis 'square'
    if (row == size(controlGamma, 1))
        xlabel('dt (ms)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('mean D^2 (deg^2)');
    else
        set(gca, 'YTickLabel', {});
    end
    legend({'D^2', 'scrambled D^2'}, 'Location', 'SouthEast');

    if ((defaults.gamma == controlGamma(iVar)) && ...
            (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'Color', [1 0 0], 'FontSize', 10);
    else
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'FontSize', 10);
    end
    drawnow;
end

function plotDisplacementD1(hFig, subplotPosVectors, iVar, ...
    controlGamma, feedbakGain, defaults, timeLagsMilliseconds, ...
    displacementX, displacementY, scrambledIntervalsDisplacementX)
% Plot the linear displacement as a function of time lag
%
% Syntax:
%   plotVelocity(hFig, subplotPosVectors, iVar, controlGamma, ...
%       feedbakGain, defaults, timeAxis, velocityArcMinPerSecond)
%
% Description:
%    Plot the velocity figure.
%
% Inputs:
%    hFig                            - Figure. The figure handle.
%    subplotPosVectors               - Matrix. A matrix of structures,
%                                      wherein each structure contains
%                                      a vector for subplot positioning
%    iVar                            - Numeric. The plot index variable.
%    controlGamma                    - Matrix. The control gamma grid.
%    feedbackGain                    - Matrix. The feedback gain grid.
%    defaults                        - Struct. The defaults structure.
%    timeLagsMilliseconds            - Vector. A vector containing the
%                                      lag in milliseconds.
%    displacementX                   - Vector. A vector containing the
%                                      displacement along the X axis.
%    displacementY                   - Vector. A vector containing the
%                                      displacement along the Y axis.
%    scrambledIntervalsDisplacementX - Vector. A vector containing the
%                                      scrambled interval displacements
%                                      along the X axis.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    [row, col] = ind2sub(size(controlGamma), iVar);
    figure(hFig);
    subplot('Position', subplotPosVectors(row, col).v);

    plot(timeLagsMilliseconds, displacementX, 'r-', 'LineWidth', 2);
    hold on;
    plot(timeLagsMilliseconds, displacementY, 'b-', 'LineWidth', 2);
    plot(timeLagsMilliseconds, scrambledIntervalsDisplacementX, 'k--', ...
        'LineWidth', 1.5);
    hold off;
    grid on
    set(gca, 'XLim', [0 500], 'YLim', [0 8], 'XTick', [0:100:1000], ...
        'YTick', [0:10], 'XScale', 'linear', 'YScale', 'linear', ...
        'FontSize', 12);
    axis 'square'
    if (row == size(controlGamma, 1))
        xlabel('dt (ms)');
    else
        set(gca, 'XTickLabel', {});
    end
    if (col == 1)
        ylabel('mean D (min arc)');
    else
        set(gca, 'YTickLabel', {});
    end

    if ((defaults.gamma == controlGamma(iVar)) && ...
            (defaults.feedback == feedbakGain(iVar)))
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'Color', [1 0 0], 'FontSize', 10);
    else
        title(sprintf('control:%0.2f, feedback:%0.2f', controlGamma(iVar), ...
            feedbakGain(iVar)), 'FontSize', 10);
    end
    drawnow;
end
