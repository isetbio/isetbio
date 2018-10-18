function t_simplePhotocurrentComputation
% Illustrate how to compute photocurrent responses to a dynamic stimulus.
%
% Syntax:
%   t_simplePhotocurrentComputation
%
% Description:
%    Simple script that demonstrates how to compute photocurrent responses
%    to a dynamic stimulus
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

% History:
%    xx/xx/18  NPC  ISETBIO Team, 2018
%    06/30/18  npc  Wrote it.
%    10/18/18  JNM  Formatting

% Set up scene
sparams.fov = 2;          % scene field of view in degrees
sparams.distance = 0.57;  % viewing distance in meters

% Gabor parameters
P(1:2) = harmonicP;
P(1).contrast = 0;

% Stimulus temporal properties
stimDurationSeconds = 200 / 1000;  % 200 milliseconds
stimRefreshRate = 10 / 1000;       % 100 Hz display
stimSampleTimes = (0:ceil(stimDurationSeconds / stimRefreshRate)) ...
    * stimRefreshRate;
stimModulationEnvelope = ieScale(fspecial('gaussian', ...
    [1, numel(stimSampleTimes)], 3), 0, 1);

% Make an optical image sequence of the stimulus
ois = oisCreate('harmonic', 'blend', stimModulationEnvelope, ...
    'sampleTimes', stimSampleTimes, ...
    'testParameters', P, ...
    'sceneParameters', sparams);

% Create cone mosaic
cMosaic = coneMosaic('center', [0, 0], 'whichEye', 'left');

% Set the field of view (degrees)
cMosaic.setSizeToFOV(sparams.fov);

% Add photon noise
cMosaic.noiseFlag = 'random';

% Set outer segment to be computed with linear filters
cMosaic.os = osLinear;

% Set the desired integration time, here 5 milliseconds
cMosaic.integrationTime = 5 / 1000;

% Compute the # of eye movements
% Eye movements are sampled according to the integration time of the cone
% mosaic. So if, the integration time is 5 msec, and the OIsequence is 1
% second long, we need to generate 200 eye movements. To compute this use
% the maxEyeMovementsNumGivenIntegrationTime() method of the oiSequence.
eyeMovementsNum = ...
    ois.maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Generate eye eyemovements for nTrials
nTrials = 30;
emPaths = cMosaic.emGenSequence(eyeMovementsNum, ...
    'nTrials', nTrials, ...
    'microsaccadetype', 'none');

% Compute absorptions and photocurrents
[absorptions, photocurrents, interpFilters, meanCur] = ...
    cMosaic.computeForOISequence(ois, ...
        'currentFlag', true, 'emPaths', emPaths);

% Visualize the mean response of the best responding L-, M- and S-cones.
% Compute absorptions & current per cone class as cones x tPoints x trials.
absorptions1 = permute(absorptions, [2 3 4 1]);
current1 = permute(photocurrents, [2 3 4 1]);
[nRows, mCols, tBins, nTrials] = size(absorptions1);
absorptions1 = reshape(absorptions1, [nRows * mCols, tBins, nTrials]);
current1 = reshape(current1, [nRows * mCols, tBins, nTrials]);

% Find indices of cones for each of the 3 cone types
lcones = cMosaic.pattern == 2;
mcones = cMosaic.pattern == 3;
scones = cMosaic.pattern == 4;

% Find indices of the best responding L, M ans S-cone
l_absorptions = absorptions1(lcones, :, :);
m_absorptions = absorptions1(mcones, :, :);
s_absorptions = absorptions1(scones, :, :);
l_curr = current1(lcones, :, :);
m_curr = current1(mcones, :, :);
s_curr = current1(scones, :, :);

[~, maxIndex] = max(l_absorptions(:));
[maxLconeIndex, ~, ~] = ind2sub(size(l_absorptions), maxIndex);
[~, maxIndex] = max(m_absorptions(:));
[maxMconeIndex, ~, ~] = ind2sub(size(m_absorptions), maxIndex);
[~, maxIndex] = max(s_absorptions(:));
[maxSconeIndex, ~, ~] = ind2sub(size(s_absorptions), maxIndex);

% Compute mean over reps of best responding L-cone
meanLconeResponse = squeeze(mean(l_absorptions(maxLconeIndex, :, :), 3));
stdLconeResponse = squeeze(std(l_absorptions(maxLconeIndex, :, :), 0, 3));

% Compute mean over reps of best responding M-cone
meanMconeResponse = squeeze(mean(m_absorptions(maxMconeIndex, :, :), 3));
stdMconeResponse = squeeze(std(m_absorptions(maxMconeIndex, :, :), 0, 3));

% Compute mean over reps of best responding S-cone
meanSconeResponse = squeeze(mean(s_absorptions(maxSconeIndex, :, :), 3));
stdSconeResponse = squeeze(std(s_absorptions(maxSconeIndex, :, :), 0, 3));

% Obtain the response time axes in milliseconds
responseTimeAxis = cMosaic.timeAxis * 1000;
tTicks = (0:stimRefreshRate:10) * 1000;

% Render results
hFig = figure(10);
clf;
set(hFig, 'Position', [10 10 750 1100]);
subplot(3, 1, 1);
% Check ois time axis, i.e. number of frames
stairs(ois.timeAxis * 1000, ois.modulationFunction, 'k', 'LineWidth', 1.5)
set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], ...
    'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
ylabel('stimulus modulation', 'FontWeight', 'bold')
xlabel('time (ms)', 'FontWeight', 'bold')

subplot(3, 1, 2)
plot(responseTimeAxis, meanLconeResponse , 'r-', 'LineWidth', 1.5);
hold on;
plot(responseTimeAxis, meanLconeResponse+stdLconeResponse , 'r:', ...
    'LineWidth', 1.0);
plot(responseTimeAxis, meanLconeResponse-stdLconeResponse , 'r:', ...
    'LineWidth', 1.0);

plot(responseTimeAxis, meanMconeResponse, 'g-', 'LineWidth', 1.5);
hold on;
plot(responseTimeAxis, meanMconeResponse+stdMconeResponse , 'g:', ...
    'LineWidth', 1.0);
plot(responseTimeAxis, meanMconeResponse-stdMconeResponse , 'g:', ...
    'LineWidth', 1.0);

plot(responseTimeAxis, meanSconeResponse, 'b-', 'LineWidth', 1.5);
hold on;
plot(responseTimeAxis, meanSconeResponse+stdSconeResponse , 'b:', ...
    'LineWidth', 1.0);
plot(responseTimeAxis, meanSconeResponse-stdSconeResponse , 'b:', ...
    'LineWidth', 1.0);

set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], ...
    'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
xlabel('time (ms)', 'FontWeight', 'bold')
ylabel('isomerizations', 'FontWeight', 'bold');

subplot(3, 1, 3);
plot(responseTimeAxis, squeeze(mean(l_curr(maxLconeIndex, :, :), 3)), ...
    'r', 'LineWidth', 1.5);
hold on
plot(responseTimeAxis, squeeze(mean(m_curr(maxMconeIndex, :, :), 3)), ...
    'g', 'LineWidth', 1.5);
plot(responseTimeAxis, squeeze(mean(s_curr(maxSconeIndex, :, :), 3)), ...
    'b', 'LineWidth', 1.5);
set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]);
set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], ...
    'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
xlabel('time (ms)', 'FontWeight', 'bold')
ylabel('photocurrent (pAmps)', 'FontWeight', 'bold');

% Visualize the oiSequence as a montage
% ois.visualize('montage');

end