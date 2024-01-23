function t_fixationalEyeMovementsCharacterize
% Explore how micro-saccade strategies affect emPath of a fixationalEM.
%
% Syntax:
%   t_fixationalEyeMovementsCharacterize
%
% Description:
%    Computes key characteristics of emPaths and explores how these differ
%    for different micro-saccade strategies the @fixationalEM object.
%
%    This function can take a while to run. You can flip useparfor to true
%    and speed up the calculation - if your machine is set up to run in
%    that mode. To do that, your startup must execute without a response!
%
%    The examined characteristics of the emPaths are:
%       - velocity
%       - fixation span
%       - power spectal density
%       - displacement analysis
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
%    07/17/18  jnm  Formatting

    close all

    % Generate 512 3-second emPaths with a sample time of 1 msec
    emDurationSeconds = 3;
    sampleTimeSeconds = 1 / 1000;
    nTrials = 512;

    % Examine different microsaccade types
    microSaccadeTypes = {'none', 'heatmap/fixation based', 'stats based'};

    for typeIndex = 1:numel(microSaccadeTypes)
        fixEMobj = fixationalEM();
        fixEMobj.randomSeed = 1;
        fixEMobj.microSaccadeType = microSaccadeTypes{typeIndex};

        % Compute the emPaths
        computeVelocity = true;
        fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
            nTrials, computeVelocity, 'useParfor', false);

        % Compute displacement and spectral analysis
        [spectrum, displacement] = performAnalyses(fixEMobj);

        % Visualize all aspects
        visualizeAnalysis(fixEMobj, spectrum, displacement, ...
            sprintf('micro-saccade type: ''%s''', ...
            microSaccadeTypes{typeIndex}));
    end

end

function [spectrum, displacement] = performAnalyses(fixEMobj)
% Analyze the emPaths
%
% Syntax:
%   [spectrum, displacement] = performAnalyses(fixEMobj)
%
% Description:
%    Analyze the emPaths
%
% Inputs:
%    fixEMobj     - Object. A fixational eye movements object.
%
% Outputs:
%    spectrum     - Struct. A structure containing such data as the vectors
%                   for frequency axis, the X power spectral density, and
%                   the Y power spectral density.
%    displacement - Struct. A structure that contains a number of vectors
%                   that deal with displacement.
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
        [dispX, ~, scrambledIntervalsDispX, ~, timeLagsMilliseconds] = ...
            fixEMobj.performDisplacementAnalysis(xPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D1');
        [dispY, ~, scrambledIntervalsDispY, ~, ~] = ...
            fixEMobj.performDisplacementAnalysis(yPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D1');

        [dispX2, ~, scrambledIntervalsDispX2, ~, ~] = ...
            fixEMobj.performDisplacementAnalysis(xPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D2');
        [dispY2, ~, scrambledIntervalsDispY2, ~, ~] = ...
            fixEMobj.performDisplacementAnalysis(yPosArcMin', ...
            fixEMobj.timeAxis, 'mode', 'D2');

        % Perform spectral analysis
        samplingRate = 1 / timeStepSeconds;
        [pSDX, frequencyAxis] = ...
            periodogram(xPosArcMin, [], [], samplingRate);
        [pSDY, ~] = periodogram(yPosArcMin, [], [], samplingRate);

        % Preallocate memory for all trials
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
        end

        % Enter trial data
        displacementX(trialNo, :) = dispX;
        displacementY(trialNo, :) = dispY;
        scrambledIntervalsDisplacementX(trialNo, :) = ...
            scrambledIntervalsDispX;
        scrambledIntervalsDisplacementY(trialNo, :) = ...
            scrambledIntervalsDispY;

        displacementX2(trialNo, :) = dispX2;
        displacementY2(trialNo, :) = dispY2;
        scrambledIntervalsDisplacementX2(trialNo, :) = ...
            scrambledIntervalsDispX2;
        scrambledIntervalsDisplacementY2(trialNo, :) = ...
            scrambledIntervalsDispY2;

        powerSpectralDensityX(trialNo, :) = pSDX;
        powerSpectralDensityY(trialNo, :) = pSDY;
    end % trialNo

    % Mean over trials spectrum
    spectrum.frequencyAxis = frequencyAxis;
    spectrum.powerSpectralDensityX = ...
        squeeze(mean(powerSpectralDensityX, 1));
    spectrum.powerSpectralDensityY = ...
        squeeze(mean(powerSpectralDensityY, 1));

    % Mean over all trials, displacement results
    displacement.timeLagsMilliseconds = timeLagsMilliseconds;
    displacement.x = squeeze(mean(displacementX, 1));
    displacement.y = squeeze(mean(displacementY, 1));
    %displacement = sqrt(displacementX .^ 2 + displacementY .^ 2);
    displacement.xScrambledIntervals = ...
        squeeze(mean(scrambledIntervalsDisplacementX, 1));
    displacement.yScrambledIntervals = ...
        squeeze(mean(scrambledIntervalsDisplacementY, 1));

    displacementX2 = squeeze(mean(displacementX2, 1));
    displacementY2 = squeeze(mean(displacementY2, 1));
    scrambledIntervalsDisplacementX2 = ...
        squeeze(mean(scrambledIntervalsDisplacementX2, 1));
    scrambledIntervalsDisplacementY2 = ...
        squeeze(mean(scrambledIntervalsDisplacementY2, 1));

    displacement.xD2degs = displacementX2 / (60 ^ 2);
    displacement.yD2degs = displacementY2 / (60 ^ 2);
    displacement.D2degs = displacement.xD2degs + displacement.yD2degs;

    displacement.xD2degsScrambedIntervals = ...
        scrambledIntervalsDisplacementX2 / (60 ^ 2);
    displacement.yD2degsScrambedIntervals = ...
        scrambledIntervalsDisplacementY2 / (60 ^ 2);
    displacement.D2degsScrambedIntervals = ...
        displacement.xD2degsScrambedIntervals + ...
        displacement.yD2degsScrambedIntervals;
end

function visualizeAnalysis(fixEMobj, spectrum, displacement, figureName)
% Display the analysis
%
% Syntax:
%   visualizeAnalysis(fixEMobj, spectrum, displacement, figureName)
%
% Description:
%    Display the emPaths.
%
% Inputs:
%    fixEMobj     - Object. A fixational eye movement object.
%    spectrum     - Struct. A structure containing such data as the vectors
%                   for frequency axis, the X power spectral density, and
%                   the Y power spectral density.
%    displacement - Struct. A structure that contains a number of vectors
%                   that deal with displacement.
%    figureName   -
%
% Outputs:
%    hFig         - Figure. The figure containing the visualization.
%
% Optional key/value pairs:
%    None.
%
    emPosRange = [-25 25];
    emPosDelta = 1;

    nTrials = size(fixEMobj.emPosArcMin, 1);
    
    % Plot the results
    hFig = figure();
    clf;
    set(hFig, 'Name', figureName, 'Position', [0 0 1150 950], ...
        'Color', [1 1 1]);

    % emPlots - position (last trial)
    subplot(3, 3, [1 2]);
    xPosArcMin = squeeze(fixEMobj.emPosArcMin(end, :, 1));
    yPosArcMin = squeeze(fixEMobj.emPosArcMin(end, :, 2));
    plot(fixEMobj.timeAxis, xPosArcMin, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(fixEMobj.timeAxis, yPosArcMin, 'b-', 'LineWidth', 1.5);
    plot(fixEMobj.timeAxis, yPosArcMin * 0, 'k-');
    hold off;
    grid on;
    timeTicks = 0:0.25:100;
    set(gca, 'XLim', [fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        'XTick', timeTicks, 'XTickLabel', {}, ...
        'YLim', emPosRange, 'FontSize', 14);
    ylabel('position (arc min)');
    legend({'x-pos', 'y-pos'})
    title('x/y eye position (last trial)');

    % Velocity (last trial)
    subplot(3, 3, [4 5]);
    plot(fixEMobj.timeAxis, squeeze(fixEMobj.velocityArcMin(end, :)), ...
        'k-', 'LineWidth', 1.5);
    hold on;
    plot(fixEMobj.timeAxis, squeeze(mean(fixEMobj.velocityArcMin, 1)), ...
        'k-', 'Color', [0.9 0.8 0.2], 'LineWidth', 3.0);
    subjectDR.mean = 30;
    subjectDR.std = 17;
    plot([fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        (subjectDR.mean - subjectDR.std) * [1 1], 'r--', 'LineWidth', 1.0);
    subjectDG.mean = 89;
    subjectDG.std = 63;
    plot([fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        (subjectDG.mean + subjectDG.std) * [1 1], 'b--', 'LineWidth', 1.0);
    hold off;
    grid on
    legend({'model (one iteration)', ...
        'model (mean across all iterations)', ...
        'Cherici et al: subjectDR', 'Cherici et al: subjectDG'}, ...
        'Location', 'NorthEast')
    set(gca, 'XLim', [fixEMobj.timeAxis(1) fixEMobj.timeAxis(end)], ...
        'XTick', timeTicks, ...
        'XTickLabel', sprintf('%2.2f\n', timeTicks), ...
        'YLim', [0 60 * 10], ...
        'YTick', [0:100:60 * 100], ...
        'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('velocity (arc min / sec)');
    title(sprintf('eye velocity'));

    subplot(3, 3, 3);
    % Last trial x/y emPath
    plot(xPosArcMin, yPosArcMin, 'k-', 'LineWidth', 1.0);
    hold on
    plot([0 0], emPosRange, 'k-');
    plot(emPosRange, [0 0], 'k-');
    hold off;
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, ...
        'XTick', [-100:10:100], 'YTick', [-100:10:100], ...
        'XTickLabel', {}, 'FontSize', 14);
    grid on;
    box on
    axis 'square'
    ylabel('y-position (arc min)');
    title('eye movement path (last trial)');

    subplot(3, 3, 6);
    hold on
    % fixation span
    cmap = brewermap(1024, 'Greys');
    maxDurationSecondsForFixationSpan = 1.0;
    [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = ...
        fixEMobj.computeFixationMap(fixEMobj.timeAxis, ...
        fixEMobj.emPosArcMin, emPosRange, emPosDelta, ...
        'maxDurationSeconds', maxDurationSecondsForFixationSpan);

    contourf(fixationMapSupportX, fixationMapSupportY, ...
        fixationMap, 0:0.05:1, 'LineColor', [.5 0.5 0.5]);
    hold on;
    plot([0 0], emPosRange, 'k-');
    plot(emPosRange, [0 0], 'k-');
    plot(fixationMapSupportX, ...
        emPosRange(1) + fixationMapXSlice * emPosRange(2) * 0.9, '-', ...
        'Color', [1 0 0], 'LineWidth', 1.5);
    plot(emPosRange(1) + fixationMapYSlice * emPosRange(2) * 0.9, ...
        fixationMapSupportY, '-', 'Color', [0 0 1], 'LineWidth', 1.5);

    hold off;
    set(gca, 'YLim', emPosRange, 'XLim', emPosRange, ...
        'XTick', [-100:10:100], 'YTick', [-100:10:100], 'FontSize', 14);
    colormap(cmap);
    % hh = colorbar('North', 'Ticks', [0:0.2:1.0], ...
    %    'TickLabels', {'.0', '.2', '.4', '.6', '.8', '1.'}, ...
    %    'Color', [.4 .4 .4]);
    xlabel('x-position (arc min)');
    ylabel('y-position (arc min)');
    grid on;
    box on;
    axis 'square';
    axis 'xy'
    title(sprintf(['fixation span during\ninitial %1.1f seconds ' ...
        '(%d trials)'], maxDurationSecondsForFixationSpan, nTrials));

    % Spectral density plots
    subplot(3, 3, 7);
    plot(spectrum.frequencyAxis, 10 * ...
        log10(spectrum.powerSpectralDensityX), 'r-', 'LineWidth', 1.5);
    hold on
    plot(spectrum.frequencyAxis, 10 * ...
        log10(spectrum.powerSpectralDensityY), 'b-', 'LineWidth', 1.5);
    hold off
    set(gca, 'XScale', 'log', 'XLim', [1 300], 'YLim', [-50 10], ...
        'XTick', [1 3 10 30 100 300], 'YTick', [-100:10:100], ...
        'FontSize', 14);
    grid on;
    legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest')
    axis 'square'
    xlabel('frequency (Hz)')
    ylabel('power density (dB)')
    title(sprintf('power spectral density\n(%d trials)', nTrials));

    % Displacement^2 plot
    subplot(3, 3, 8);
    plot(displacement.timeLagsMilliseconds, displacement.D2degs, ...
        'k-', 'LineWidth', 1.5);
    hold on;
    plot(displacement.timeLagsMilliseconds, ...
    displacement.D2degsScrambedIntervals, 'k--', 'LineWidth', 1.5);
hold off;
    grid on
    minD2 = min([min(displacement.D2degs(:)), ...
        min(displacement.D2degsScrambedIntervals(:))]);
    D2range = [minD2 minD2 * 1000];
    set(gca, 'XLim', [1 1000], 'YLim', D2range, ...
        'XTick', [1 3 10 30 100 300 1000], ...
        'YTick', [1e-4 1e-3 1e-2 1e-1], ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    axis 'square'
    xlabel('time interval (ms)')
    ylabel('mean displacement^2 (deg^2)');
    legend({'emPath', 'scrambled'}, 'Location', 'SouthEast')
    title(sprintf('mean squared displacement\n(%d trials)', nTrials));

    % Linear displacement plot
    subplot(3, 3, 9);
    plot(displacement.timeLagsMilliseconds, displacement.x, 'r-', ...
        'LineWidth', 2);
    hold on;
    plot(displacement.timeLagsMilliseconds, displacement.y, 'b-', ...
        'LineWidth', 2);
    plot(displacement.timeLagsMilliseconds, ...
        displacement.xScrambledIntervals, 'k--', 'LineWidth', 1.5);
    hold off;
    grid on
    minD = min([min(displacement.x(:)), ...
        min(displacement.xScrambledIntervals(:))]);
    Drange = [minD minD * 100];

    set(gca, 'XLim', [2 1000], 'YLim', Drange, ...
        'XTick', [3 10 30 100 300 1000], ...
        'YTick', [0.03 0.1 0.3 1 3 10], ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    axis 'square'
    xlabel('time interval (ms)')
    ylabel('mean displacement (min arc)')
    legend({'emPath (x)', 'emPath (y)', 'scrambled (x)'}, ...
        'Location', 'SouthEast');
    title(sprintf('mean displacement\n(%d trials)', nTrials));
    drawnow;
end
