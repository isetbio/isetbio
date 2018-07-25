function t_fixationalEyeMovementsToIsomerizations(varargin)
% Show eye movements on cone mosaic and sequence of cone activations
%
% Syntax:
%   t_fixationalEyeMovementsToIsomerizations()
%
% Description:
%    Shows how to integrate fixational eye movements with a cone mosaic and
%    an oiSequence, and illustrates the different time bases for the
%    stimulus modulation, eye-movements and isomerizations.
%
%    This tutorial simulates a small uniform square stimulus (0.05 deg)
%    whose luminance is modulated over time, and which is transported over
%    the cone mosaic due to eye movements. The tutorial produces two
%    figures: The first one shows the cone mosaic with one eye movement
%    path superimposed on top of it.
%
%    The second figure is a multi-panel figure and displays the following:
%       - the stimulus spatial profile (top left)
%       - 2D mosaic activations (isomerizations) at different epochs during
%         the course of the trial
%       - Space (cone) x Time plot of the isomerization signals. showing
%         how the mosaic response evolves during the entire trial
%       - Time plot of the x- and y- eye position during the entire trial.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    'fov degrees'        - Numeric. The number of degrees in the field of
%                           view. Default is 0.25 degrees
%    'integration time'   - Numeric. The integration time in seconds.
%                           Default is 1ms (1e-3 sec).
%    'resampling factor'  - Numeric. The resampling factor. Default is 15
%    'use parfor'         - Boolean. Whether to compute using a parfor loop
%                           (to take advantage of multiple processors) of
%                           to use a single processor. Default false, i.e.,
%                           use a single processor.
%

% History
%    02/06/18  npc  Wrote it.
%    02/07/18  npc  Comments.
%    07/11/18  jnm  Formatting. Written as a function, so formatted as one.

% Examples:
%{
    t_fixationalEyeMovementsToIsomerizations;
%}

%%
p = inputParser;
varargin = ieParamFormat(varargin);

% Cone hexagonal mosaic params. Here we specify that the mosaic will have a
% field of view of 0.25 degrees and integration time of 1ms The resampling
% factor of 15  means that the hex locations will be computed by sampling a
% rectangular grid whose pixel size is 1/15th of the pixel size of the
% default rectangular mosaic. The higher the resamplingFactor, the closer
% the cone positions wil be to the nodes of a perfectly hexagonal grid.

p.addParameter('fovdegs', 0.25, @isscalar);
p.addParameter('integrationtime', 1e-3, @isscalar);
p.addParameter('resamplingfactor', 15, @isscalar);
p.addParameter('useparfor', false, @islogical);

p.parse(varargin{:});
fovDegs = p.Results.fovdegs;
integrationTime = p.Results.integrationtime;
resamplingFactor = p.Results.resamplingfactor;
useParfor = p.Results.useparfor;

%% Stimulus params.
% Here we are specifying a square stimulus that is 0.05 degs in width, has
% a max luminance of 100 cd/m2 and 0.5 second duration. The stimulus
% refresh interval is 20 milliseconds, and the temporal modulation has a
% pulse envelope with an amplitude of 1.0 between frames 5 and 20 (i.e.,
% 100-400 msecs) and 0.1 for the remaining time.
sceneParams = struct('fov', 0.05, 'luminance', 100);
stimDurationSeconds = 0.5;
stimRefreshInterval = 20 / 1000;
stimFrames = round(stimDurationSeconds / stimRefreshInterval);
stim.modulation = 0.1 + zeros(1, stimFrames);
stim.modulation(5:20) = 1.0;
stim.timeAxis = stimRefreshInterval * ((1:length(stim.modulation)) - 1);

% Generate an oiSequence for the pulse stimulus. Here we generate a
% sequence of optical images in which each optical image corresponds to a
% different frame of the stimulus sequence.
theOIsequence = oisCreate('impulse', 'add', stim.modulation, ...
    'sampleTimes', stim.timeAxis, 'sceneParameters', sceneParams);

%% A fixational eye movement superimposed on the cone mosaic
% Generating fixational eye movements that include drift and microsaccades.
fixEMobj = fixationalEM();

% Instantiate a hexagonal cone mosaic.
cm = coneMosaicHex(resamplingFactor, 'fovDegs', fovDegs);
cm.integrationTime = integrationTime;

% Compute the number of eye movements for the given oiSequence and the
% integration time (which is also the timeSample for the eye movements)
% of the cone mosaic which will be used to compute responses.
eyeMovementsPerTrial = ...
    theOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
    cm.integrationTime);

% Compute emPaths for this mosaic for nTrials of the same oiSequence. Here
% we are fixing the random seed so as to reproduce identical eye movements
% whenever this script is run.
nTrials = 2;
tic
fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
    'nTrials', nTrials, ...
    'rSeed', 857, ...
    'useParfor', useParfor);
toc
visualizedTrial = 1;

emPlot(fixEMobj, 'emMosaic', 'cone mosaic', cm, ...
    'visualized trial', visualizedTrial);
%% Compute the mosaic response for all trials.

%  Here we use the fixational eye movements in units of cone mosaic pattern
%  size ('emPos'), which is what the  @coneMosaic.computeForOISequence()
%  method expects.
fprintf(['Compute mosaic cone responses (isomerizations) to the ' ...
    'background sequence, taking eye movements into account.\n']);
[isomerizations, ~, ~, ~] = ...
    cm.computeForOISequence(theOIsequence, ...
    'emPaths', fixEMobj.emPos, ...
    'seed', 1, ...
    'interpFilters', [], ...
    'meanCur', [], ...
    'currentFlag', false);

% The returned isomerizations are in photons/integration time.
% Convert to isomerization rate by dividing with the integrationTime
isomerizationRate = isomerizations / cm.integrationTime;

% Form a struct with the response data to visualize
response.timeAxis = cm.timeAxis();
response.isomerizationRate = ...
    squeeze(isomerizationRate(visualizedTrial, :, :));
response.isomerizationRange = [0 max(response.isomerizationRate(:))];

% Extract the time range for all the signals
timeLimits = [0, max([max(stim.timeAxis), ...
    max(response.timeAxis), max(fixEMobj.timeAxis)])];

hFig = vcNewGraphWin;
set(hFig, 'Position', [0 0 1.0 0.91]);

%% Plot the stimulus modulation function and superimpose the stimulus
% frames (gray lines)
subplot(3, 6, 1);
hold on;
% Superimpose stimulus frames
for k = 1:numel(stim.timeAxis)
    plot(stim.timeAxis(k) * [1 1], stim.modulation(k) * [0 1], 'k-', ...
        'Color', [0.1 0.7 0.3], 'LineWidth', 1.5);
end
stairs(theOIsequence.timeAxis, theOIsequence.modulationFunction, 'r-', ...
    'LineWidth', 1.5);
set(gca, 'XLim', timeLimits, 'YLim', [0 1.1], 'FontSize', 12);
xlabel('time (sec)');
ylabel('stim modulation');
axis 'square'
title('stimulus modulation');

% Times during the response for visualizing the 2D mosaic acivation.
sampledTimes = [75 245 365 380 455] / 1000;

% Colormap for visualizing mosaic activation
activationColorMap = brewermap(1024, 'YlOrBr');

% Render the mosaic activation at the selected times
for k = 1:numel(sampledTimes)
    [~, responseFrame] = min(abs(response.timeAxis - sampledTimes(k)));
    ax = subplot(3, 6, 1 + k);
    activation = squeeze(response.isomerizationRate(:, responseFrame));
    cm.renderActivationMap(ax, activation, ...
        'signalRange', response.isomerizationRange, ...
        'mapType', 'modulated disks', ...
        'colorMap', activationColorMap);
    title(sprintf('time: %2.2f sec', response.timeAxis(responseFrame)));
    set(gca, 'FontSize', 12);
end

%% Render the mosaic activation as a cone(y)-time(x) plot.
ax = subplot(3, 6, [7 12]);
conesNum = size(response.isomerizationRate, 1);
imagesc(response.timeAxis, 1:conesNum, response.isomerizationRate);
axis 'xy';
hold on;
% Superimpose stimulus frames
for k = 1:numel(stim.timeAxis)
    plot(stim.timeAxis(k) * [1 1], [1 conesNum * stim.modulation(k)], ...
        'k-', 'Color', [0.1 0.7 0.3], 'LineWidth', 1.5);
end
% Superimpose times at which we display the 2D mosaic activation maps
for k = 1:numel(sampledTimes)
    plot(sampledTimes(k) * [1 1], [1 conesNum], 'k-v', 'LineWidth', 1.5);
end
set(gca, 'XLim', timeLimits, 'YLim', [1 conesNum]);
set(gca, 'CLim', response.isomerizationRange, 'FontSize', 12);
ylabel('cone id');
title('isomerization rates');
colormap(ax, colormap(ax, activationColorMap));
colorbar('Location', 'East', 'Ticks', 0:10000:35000, ...
    'TickLabels', {'0', '10K', '20K', '30K'});

%% Plot the x and y eye position during the course of the trial
% Superimpose the stimulus frames and the times at which we display the
% 2D activation maps
subplot(3, 6, [13 18]);
plot(fixEMobj.timeAxis, ...
    squeeze(fixEMobj.emPosMicrons(visualizedTrial, :, 1)), 'r-', ...
    'LineWidth', 2);
hold on
plot(fixEMobj.timeAxis, ...
    squeeze(fixEMobj.emPosMicrons(visualizedTrial, :, 2)), 'b-', ...
    'LineWidth', 2);
plot(fixEMobj.timeAxis, 0 * fixEMobj.timeAxis, 'k-', ...
    'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
% Superimpose stimulus frames
for k = 1:numel(stim.timeAxis)
    plot(stim.timeAxis(k) * [1 1], stim.modulation(k) * [-100 100], ...
        'k-', 'Color', [0.1 0.7 0.3], 'LineWidth', 1.5);
end
% Superimpose times at which we display the 2D mosaic activation maps
for k = 1:numel(sampledTimes)
    plot(sampledTimes(k) * [1 1], [-40 40], 'k-v', 'LineWidth', 1.5);
end
xlabel('time (sec)');
ylabel('eye position');
title('eye movements')
legend({'x-eye pos', 'y-eye pos'});
set(gca, 'XLim', timeLimits, 'YLim', [-40 40], ...
    'FontSize', 12, 'Color', squeeze(activationColorMap(1, :)));
end