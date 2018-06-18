function t_conesMosaicAbsorptionOISequenceEcc

%%t_conesMosaicAbsorptionOISequenceEcc  Illustrates how to compute  
%  responses to a dynamic stimulus.  
%
% Description:
%    Illustrate
%    
% 
% NPC, ISETBIO Team, 2018
%
% 06/16/18  npc  Wrote it.

% Set up scene
sparams.fov = 1;            % scene field of view in degrees
sparams.distance = 0.57;    % viewing distance in meters

% Gabor parameters
P(1:2) = harmonicP;
P(1).contrast = 0;
P(2).ph = 0;

% Stimulus temporal properties
stimDurationSeconds = 500/1000;   % 50 milliseconds
stimRefreshRate = 10/1000;       % 100 Hz display
stimSampleTimes = (0:ceil(stimDurationSeconds/stimRefreshRate))*stimRefreshRate;
stimModulationEnvelope = ieScale(fspecial('gaussian',[1,numel(stimSampleTimes)],3),0,1);

% Make an optical image sequence of the stimulus
ois = oisCreate('harmonic','blend', stimModulationEnvelope, ...
    'sampleTimes', stimSampleTimes, ...
    'testParameters', P, ...
    'sceneParameters', sparams);

% Visualize the oiSequence as a montage
%ois.visualize('montage');


% Create cone mosaic
saveMosaic = ~true;
mosaicFileName = sprintf('mosaic.mat');
if (saveMosaic)
    cMosaic  = coneMosaicHex(7, ...
            'name', 'test', ...
            'fovDegs', sparams.fov, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', ~true, ...
            'sConeMinDistanceFactor', 2.0, ... 
            'sConeFreeRadiusMicrons', 0.15*300, ...                   
            'spatialDensity', [0 6/10 3/10 1/10], ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...         
            'latticeAdjustmentDelaunayToleranceF', 0.05, ... 
            'maxGridAdjustmentIterations', 20, ...
            'marginF', [] ... 
        );
    save(mosaicFileName, 'cMosaic', '-v7.3');
else
    load(mosaicFileName, 'cMosaic');
end

% Apply ecc-based correction in cone efficiency
cMosaic.eccBasedConeQuantalEfficiency = ~true;

% Add photon noise
cMosaic.noiseFlag = 'none'; % 'random';

% Set the desired integration time, here 5 milliseconds
cMosaic.integrationTime = 5/1000;

% Compute the # of eye movements
% Eye movements are sampled according to the integration time of the cone
% mosaic. So if, the integration time is 5 msec, and the OIsequence is 1
% second long, we need to generate 200 eye movements. To compute this use
% the maxEyeMovementsNumGivenIntegrationTime() method of the oiSequence.
eyeMovementsNum = ois.maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Generate eye eyemovements for nTrials
nTrials = 32;
emPaths = cMosaic.emGenSequence(eyeMovementsNum, ...
    'nTrials', nTrials, ...
    'microsaccadetype', 'none');
emPaths = emPaths * 0;

fprintf('Computing no ecc correction\n');
% Compute absorptions and photocurrents
absorptions = ...
    cMosaic.computeForOISequence(ois, ...
        'currentFlag', false, ...
        'emPaths', emPaths);

[lConeIndices, mConeIndices,sConeIndices]  = indicesForCones(cMosaic);

% Find indices of the best responding L, M ans S-cone
iTrial = 1;
l_absorptions = squeeze(absorptions(iTrial,lConeIndices,:));
m_absorptions = squeeze(absorptions(iTrial,mConeIndices,:));
s_absorptions = squeeze(absorptions(iTrial,sConeIndices,:));

[~,maxIndex] = max(l_absorptions(:));
[maxLconeIndex, timeBinOfMaxResponse] = ind2sub(size(l_absorptions), maxIndex);
[~,maxIndex] = max(m_absorptions(:));
[maxMconeIndex, ~] = ind2sub(size(m_absorptions), maxIndex);
[~,maxIndex] = max(s_absorptions(:));
[maxSconeIndex, ~] = ind2sub(size(s_absorptions), maxIndex);


% Obtain the response time axes in milliseconds
responseTimeAxis = cMosaic.timeAxis * 1000;
tTicks = (0:stimRefreshRate:10)*1000;

% Render results
hFig = figure(10); clf;
set(hFig, 'Position', [10 10 750 1100]);
subplot(2,3,1);
% Check ois time axis, i.e. number of frames
stairs(ois.timeAxis*1000,ois.modulationFunction,'k', 'LineWidth', 1.5)
set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], 'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
ylabel('stimulus modulation', 'FontWeight', 'bold')
xlabel('time (ms)', 'FontWeight', 'bold')

subplot(2,3,2)
plot(responseTimeAxis,squeeze(l_absorptions(maxLconeIndex,:)) ,'r-', 'LineWidth', 1.5); hold on;
plot(responseTimeAxis,squeeze(m_absorptions(maxMconeIndex,:)),'g-', 'LineWidth', 1.5);
plot(responseTimeAxis,squeeze(s_absorptions(maxSconeIndex,:)),'b-', 'LineWidth', 1.5);

set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], 'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
xlabel('time (ms)', 'FontWeight', 'bold')
ylabel('isomerizations', 'FontWeight', 'bold');

subplot(2,3,3)
activation = squeeze(absorptions(iTrial,:,timeBinOfMaxResponse));
activation2D = cMosaic.reshapeHex1DmapToHex2Dmap(activation);
signalRange = [min(activation(:)) max(activation(:))];
cMosaic.renderActivationMap(gca, activation2D, ...
         'signalRange', signalRange, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0]);
title('Without ecc-based corrections');

% Computing again
fprintf('Computing with ecc on now\n');
cMosaic.eccBasedConeQuantalEfficiency = true;

% Compute absorptions and photocurrents
absorptions = ...
    cMosaic.computeForOISequence(ois, ...
        'currentFlag', false, ...
        'emPaths', emPaths);
    
[lConeIndices, mConeIndices,sConeIndices]  = indicesForCones(cMosaic);

% Find indices of the best responding L, M ans S-cone
iTrial = 1;
l_absorptions = squeeze(absorptions(iTrial,lConeIndices,:));
m_absorptions = squeeze(absorptions(iTrial,mConeIndices,:));
s_absorptions = squeeze(absorptions(iTrial,sConeIndices,:));

[~,maxIndex] = max(l_absorptions(:));
[maxLconeIndex, timeBinOfMaxResponse] = ind2sub(size(l_absorptions), maxIndex);
[~,maxIndex] = max(m_absorptions(:));
[maxMconeIndex, ~, ] = ind2sub(size(m_absorptions), maxIndex);
[~,maxIndex] = max(s_absorptions(:));
[maxSconeIndex, ~, ] = ind2sub(size(s_absorptions), maxIndex);


% Obtain the response time axes in milliseconds
responseTimeAxis = cMosaic.timeAxis * 1000;
tTicks = (0:stimRefreshRate:10)*1000;

% Render results

subplot(2,3,4);
% Check ois time axis, i.e. number of frames
stairs(ois.timeAxis*1000,ois.modulationFunction,'k', 'LineWidth', 1.5)
set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], 'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
ylabel('stimulus modulation', 'FontWeight', 'bold')
xlabel('time (ms)', 'FontWeight', 'bold')

subplot(2,3,5)
plot(responseTimeAxis,squeeze(l_absorptions(maxLconeIndex,:)) ,'r-', 'LineWidth', 1.5); hold on;
plot(responseTimeAxis,squeeze(m_absorptions(maxMconeIndex,:)),'g-', 'LineWidth', 1.5);
plot(responseTimeAxis,squeeze(s_absorptions(maxSconeIndex,:)),'b-', 'LineWidth', 1.5);

set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)], 'XTick', tTicks);
grid on
set(gca, 'FontSize', 12);
xlabel('time (ms)', 'FontWeight', 'bold')
ylabel('isomerizations', 'FontWeight', 'bold');

subplot(2,3,6)
activation = squeeze(absorptions(iTrial,:,timeBinOfMaxResponse));
activation2D = cMosaic.reshapeHex1DmapToHex2Dmap(activation);
cMosaic.renderActivationMap(gca, activation2D, ...
         'signalRange', signalRange, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0]);
title('With ecc-based corrections');


end