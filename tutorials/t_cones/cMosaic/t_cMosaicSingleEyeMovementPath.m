% Demo basic usage of the new @cMosaic object with eye movements
%
% Description:
%    Shows basic usage of the new cone mosaic class, @cMosaic with eye
%    movements. Here, we generate an on-axis (zero eccentricity) cMosaic object and
%    compute N noisy response instances to a static stimulus
%    under a single fixational eye movement path.
%
% See Also:
%   t_cMosaicMultipleEyeMovementPaths
%   t_cMosaicBasic
%   t_cMosaicGenerate
%   t_cMosaicEccDependentAbsorptionEfficacy
%   t_cMosaicBenchMark
%   t_cMosaicOffAxis
%   t_cMosaicFromConeMosaicHex

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate the ring rays stimulus
parms.freq = 4;
parms.contrast = 1;
parms.ph = 0;
parms.ang = 0;
parms.row = 128;
parms.col = 128;
parms.GaborFlag = 0;
scene = sceneCreate('harmonic', parms);
scene = sceneSet(scene, 'fov', 0.25);

%% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

%% Generate the mosaic
cm = cMosaic(...
    'sizeDegs', [0.25 0.25], ...      % SIZE: 1.0 degs (x) 0.5 degs (y)
    'eccentricityDegs', [0 0], ...  % ECC: (0,0)
    'integrationTime', 5/1000 ...   % integration time: 5 msec
    );

%% Generate 1 eye movement path lasting for 100 msec
eyeMovementDurationSeconds = 100/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'nTrials', 1, ...
        'randomSeed', 10);
    
%% Compute 128 noisy response instances of cone excitation response to the same eye movement path
instancesNum = 128;
[noiseFreeExcitationResponseInstances, noisyExcitationResponseInstances, ~,~,timeAxis] = cm.compute(oi, ...
        'withFixationalEyeMovements', true, ...
        'nTrials', instancesNum);

    
    
%% Visualize the mosaic
cm.visualize();

%% Visualize time-series response of a singe cone
% Lets plot responses for the cone with max noise-free response
[~,idx] = max(noiseFreeExcitationResponseInstances(:));
[~,~,targetConeID] = ind2sub(size(noiseFreeExcitationResponseInstances), idx);

figure(2); clf;
% Plot the time series response for individual instances
plot(timeAxis, squeeze(noisyExcitationResponseInstances(:,:,targetConeID)), 'k.');
hold on;
% Plot the time series response for the mean of the individual instances
plot(timeAxis, squeeze(mean(noisyExcitationResponseInstances(:,:,targetConeID),1)), 'g-', 'LineWidth', 2.0);
% Plot the noise-free time series response in red
plot(timeAxis, squeeze(noiseFreeExcitationResponseInstances(:,:,targetConeID)), 'r', 'LineWidth', 1.5);
xlabel('time (seconds)');
ylabel('excitations per integration time');
set(gca, 'FontSize', 16);



%% Visualize time series cone mosaic response (noise-free)
% Extract the eye movement path to visualize
emPathsDegs = cm.fixEMobj.emPosArcMin/60;

hFig = figure(3); clf;
set(hFig, 'Position', [100 300 1120 1060]);
activationRange = prctile(noiseFreeExcitationResponseInstances(:), [1 99]);

subplotPos = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 4, ...
       'colsNum', 5, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.03);

for timePoint = 1:size(noiseFreeExcitationResponseInstances,2)
    if (timePoint <= 20)
        row = floor((timePoint-1)/5)+1;
        col = mod(timePoint-1,5)+1;
        ax = subplot('Position', subplotPos(row,col).v);
        cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', noiseFreeExcitationResponseInstances(1,timePoint,:), ...
             'activationRange', activationRange, ...
             'crossHairsOnOpticalImageCenter', true, ...
             'currentEMposition', squeeze(emPathsDegs(1,timePoint,:)), ...
             'displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:timePoint), ...
             'noXLabel', true, 'noYLabel', true, ...
             'plotTitle',  sprintf('%2.0f msec', timeAxis(timePoint)*1000));
    end
    
end

