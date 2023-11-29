function t_fixationalEMTimeSampling
% Examine results of time sampling on resulting eye movement paths.
%
% Syntax:
%   t_fixationalEMTimeSampling
%
% Description:
%    Examine the effects of different time sampling (cone integration time)
%    on the resulting eye movement paths. Shows that the path doesn't
%    change drastically when we compute with a longer time interval up to  
%    five milliseconds. 
%    We should compare the frequency component as well. Some high frequency
%    content will be lost at longering intergation times
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
%    07/18/18  jnm  Formatting

    close all;

    % Cone mosaic params
    %
    % Use just small patch of mosaic so it runs fairly quickly.
    % 1 msec is fast enough for most calculations we'd do, look at
    % effect of slowing down from there.
    fovDegs = 0.1;  
    
    % Different integration times to examine
    integrationTimes = [1.0 2.5 4.0 5.0 10]/1000;
    
    % Stimulus params for a pulse stimulus
    sceneParams = struct('fov', fovDegs, 'luminance', 100);
    stimRefreshInterval = 50 / 1000;
    stimWeights = zeros(1,50);
    stimWeights(4) = 1;
    sampleTimes = stimRefreshInterval * ((1:length(stimWeights)) - 1);

    % Generate an oiSequence for the pulse stimulus
    theOIsequence = oisCreate('impulse','add', stimWeights, ...
        'sampleTimes', sampleTimes, 'sceneParameters', sceneParams);

    % Instantiate a hex mosaic
    resamplingFactor = 9;
    cm = coneMosaicHex(resamplingFactor, 'fovDegs', fovDegs);

    % Instantiate a fixational eye movement object
    fixEMobj = fixationalEM();
    
    % Generate 2 trials to make it run faster
    nTrials = 2;

    % Set up figure
    vcNewGraphWin([], 'wide');
    colors = [1 0 0; 0 0 1; 0 0 0; 0.4 0.4 0.4; 0 0 0];
    legends = {};

    for intTimeIndex = 1:numel(integrationTimes)
        % Set the mosaic's integration time
        % (this is also the time sample for the emPaths)
        cm.integrationTime = integrationTimes(intTimeIndex);
        legends{numel(legends) + 1} = sprintf('int. time: %2.0f ms', ...
            round(integrationTimes(intTimeIndex) * 1000));

        % Compute the number of eye movements given this integration time
        % and oiSequence
        eyeMovementsPerTrial = ...
            theOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
            cm.integrationTime);

        % Compute emPath for this mosaic using same random seed, so we can
        % compare the effects of different time sampling.
        fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
            'nTrials', nTrials, 'computeVelocity', true, 'rSeed', 1);

        % Visualize the first trial emPath and velocity
        visualizedTrial = 1;
        subplot(3,1,1);
        hold on
        plot(fixEMobj.timeAxis * 1000, ...
            squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,1)), 's-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.8 0.8 0.8], ...
            'Color', squeeze(colors(intTimeIndex,:)));
        legend(legends);
        xlabel('time (ms)')
        ylabel('x-position (arc min)');
        grid on
        set(gca, 'FontSize', 14);

        subplot(3,1,2);
        hold on
        plot(fixEMobj.timeAxis * 1000, ...
            squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,2)), 's-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.8 0.8 0.8], ...
            'Color', squeeze(colors(intTimeIndex,:)));
        legend(legends);
        xlabel('time (ms)')
        ylabel('y-position (arc min)');
        grid on
        set(gca, 'FontSize', 14);

        subplot(3,1,3);
        hold on
        plot(fixEMobj.timeAxis * 1000, ...
            squeeze(fixEMobj.velocityArcMin(visualizedTrial,:)), 's-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.8 0.8 0.8], ...
            'Color', squeeze(colors(intTimeIndex,:)));
        legend(legends);
        xlabel('time (ms)')
        ylabel('velocity (arc min / sec)');
        grid on
        set(gca, 'FontSize', 14);
    end

end
