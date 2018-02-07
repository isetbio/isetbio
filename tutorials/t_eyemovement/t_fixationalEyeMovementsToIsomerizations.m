function t_fixationalEyeMovementsToIsomerizations
% Shows how to integrate fixational eye movements with a cone mosaic and 
% an oiSequence, and illustrates the different time bases for the stimulus
% modulation, eye-movements and isomerizations.

% History
%   02/06/18  npc  Wrote it.

    close all
    
    % cone mosaic params
    fovDegs = 0.25;
    integrationTime = 1/1000;
    resamplingFactor = 15;
    
    % Stimulus params for a pulse stimulus
    stimFOVdegs = 0.05;
    sceneParams = struct('fov', stimFOVdegs, 'luminance', 100);
    stimDurationSeconds = 0.5;
    stimRefreshInterval = 20/1000;
    stimFrames = round(stimDurationSeconds/stimRefreshInterval);
    stim.modulation = 0.1+zeros(1,stimFrames); 
    stim.modulation(5:20) = 1.0;
    stim.timeAxis = stimRefreshInterval*((1:length(stim.modulation))-1);
    
    % Generate an oiSequence for the pulse stimulus
    theOIsequence = oisCreate('impulse','add', stim.modulation, ...
        'sampleTimes', stim.timeAxis, ...
        'sceneParameters', sceneParams);
    
    % Instantiate a fixational eye movement object
    fixEMobj = fixationalEM();
    nTrials = 2;
    
    % Instantiate a hex mosaic 
    cm = coneMosaicHex(resamplingFactor, 'fovDegs', fovDegs);
    cm.integrationTime = integrationTime;
    
    % Compute the number of eye movements for this integration time and oiSequence
    eyeMovementsPerTrial = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
        
    % Compute emPath for this mosaic using same random seed, so we can
    % compare the effects of different time sampling.
    fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
        'nTrials', nTrials, ...
        'computeVelocity', true, ...
        'rSeed', 857);
    
    % Visualize the emPath for one trial on top of the cone mosaic
    % Here we use the 'emPosMicrons' data
    visualizedTrial = 1;
    hFig = vcNewGraphWin;
    set(hFig, 'Position', [0 0 0.5 0.75]);
    cm.visualizeGrid(...
        'axesHandle', gca, ...
        'overlayEMpathmicrons', squeeze(fixEMobj.emPosMicrons(visualizedTrial,:,:)), ...
        'overlayNullSensors', false, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', true,...
        'generateNewFigure', false);
    xyRangeMeters = max(abs(fixEMobj.emPosMicrons(:)))*[-1 1]*1e-6;
    ticks = [-200:20:200]*1e-6;
    tickLabels = sprintf('%2.0f\n', ticks*1e6);
    set(gca, 'XLim', xyRangeMeters, 'YLim', xyRangeMeters, ...
        'XTick', ticks, 'YTick', ticks, 'XTickLabels', tickLabels, 'YTickLabels', tickLabels);
    xlabel('space (microns)');

    
    % Compute the mosaic response. Here we use the 'emPos' data, which
    % contain the eye movement paths in units of patternSampleSize
    fprintf('Compute mosaic cone responses (isomerizations) to the background sequence, taking eye movements into account.\n');
    [isomerizations, ~,~,~] = ...
    cm.computeForOISequence(theOIsequence, ...
        'emPaths', fixEMobj.emPos, ...
        'seed', 1, ...
        'interpFilters', [], ...
        'meanCur', [], ...
        'currentFlag', false);
    
    isomerizationRate = isomerizations / cm.integrationTime;
    
    % Visualize isomerizations for this trial
    response.timeAxis = cm.timeAxis();
    response.isomerizationRate = squeeze(isomerizationRate(visualizedTrial, :, :));
    response.isomerizationRange = [0 max(response.isomerizationRate(:))];

    timeLimits = [0 max([max(stim.timeAxis) max(response.timeAxis) max(fixEMobj.timeAxis)])];
    
    hFig = vcNewGraphWin;
    set(hFig, 'Position', [0 0 1.0 0.91]);
    subplot(3,6,1);
    stairs(stim.timeAxis, stim.modulation, 'r-', 'LineWidth', 1.5);
    hold on;
    for k = 1:numel(stim.timeAxis)
        plot(stim.timeAxis(k)*[1 1], stim.modulation(k)*[0 1], 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    set(gca, 'XLim', timeLimits, 'YLim', [0 1.1]);
    xlabel('time (sec)');
    ylabel('stim modulation');
    axis 'square'
    title('stimulus modulation');
    set(gca, 'FontSize', 12);
    
    sampledTimes = [70 200 360 380 450]/1000;
    for k = 1:numel(sampledTimes)
        [~,responseFrame] = min(abs(response.timeAxis-sampledTimes(k)));
        ax = subplot(3,6,1+k);
        activation = squeeze(response.isomerizationRate(:,responseFrame));
        cm.renderActivationMap(ax, activation, ...
            'signalRange', response.isomerizationRange, ...
            'colorMap', brewermap(1024, 'YlOrRd'));
        title(sprintf('time: %2.2f sec', response.timeAxis(responseFrame)));
        set(gca, 'FontSize', 12);
    end
    
    ax = subplot(3,6,[8 12]);
    conesNum = size(response.isomerizationRate,1);
    imagesc(response.timeAxis, 1:conesNum, response.isomerizationRate);
    axis 'xy';
    hold on;
    for k = 1:numel(stim.timeAxis)
        plot(stim.timeAxis(k)*[1 1], [1 conesNum*stim.modulation(k)], 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
        
    set(gca, 'XLim', timeLimits, 'YLim', [1 conesNum]);
    set(gca, 'CLim', response.isomerizationRange);
    ylabel('cone id');
    title('isomerization rates');
    set(gca, 'FontSize', 12);
    colormap(ax, colormap(ax, brewermap(1024, 'YlOrRd')));
    colorbar('Location', 'East', 'Ticks', [0:10000:35000], 'TickLabels', {'0', '10K', '20K', '30K'});
    
    
    subplot(3,6,[14 18]);
    plot(fixEMobj.timeAxis, squeeze(fixEMobj.emPosMicrons(visualizedTrial,:,1)), 'r-', 'LineWidth', 1.5); hold on
    plot(fixEMobj.timeAxis, squeeze(fixEMobj.emPosMicrons(visualizedTrial,:,2)), 'b-', 'LineWidth', 1.5);
    plot(fixEMobj.timeAxis, 0*fixEMobj.timeAxis, 'k-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
    for k = 1:numel(stim.timeAxis)
        plot(stim.timeAxis(k)*[1 1], stim.modulation(k)*[-100 100], 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    for k = 1:numel(sampledTimes)
        plot(sampledTimes(k)*[1 1], [-100 100], 'k-v', 'LineWidth', 1.5);
    end
    xlabel('time (sec)');
    title('eye movements')
    legend({'x-eye pos', 'y-eye pos'});
    set(gca, 'XLim', timeLimits, 'YLim', [-100 100]);
    set(gca, 'FontSize', 12);
    
end


