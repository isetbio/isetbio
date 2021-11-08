% Demo basic usage of the new @cMosaic object with eye movements
%
% Description:
%    Shows basic usage of the new cone mosaic class, @cMosaic with eye
%    movements. Here, we generate an on-axis (zero eccentricity) cMosaic object and
%    compute N noisy response instances to a static stimulus
%    with each instance occuring under a different fixational eye movement path.
%
% See Also:
%   t_cMosaicSingleEyeMovementPath
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

%% Generate 64 eye movement paths, each lasting for 100 msec
eyeMovementDurationSeconds = 100/1000;
instancesNum = 64;
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'nTrials', instancesNum, ...
        'randomSeed', 10);
    
%% Compute noisy response instances of cone excitation response to each of these eye movement paths
[noiseFreeExcitationResponseInstances, noisyExcitationResponseInstances, ~,~,timeAxis] = cm.compute(oi, ...
        'withFixationalEyeMovements', true);

    
    
%% Visualize the mosaic
cm.visualize();

%% Visualize time-series responses of a singe cone to each of the eye movement path
% Lets plot responses for the cone with max noise-free response
[~,idx] = max(noiseFreeExcitationResponseInstances(:));
[~,~,targetConeID] = ind2sub(size(noiseFreeExcitationResponseInstances), idx);


hFig = figure(2); clf;
set(hFig, 'Position', [100 300 1120 1060]);
subplotPos = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 8, ...
       'colsNum', 8, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.03);
   
for instanceIndex = 1:instancesNum
    row = floor((instanceIndex-1)/8)+1;
    col = mod(instanceIndex-1,8)+1;
    ax = subplot('Position', subplotPos(row,col).v);
    % Plot the noise-free response in red
    plot(timeAxis, squeeze(noiseFreeExcitationResponseInstances(instanceIndex,:,targetConeID )), 'r-', 'LineWidth', 1.5);
    hold on;
    % Plot the noisy series response in black
    plot(timeAxis, squeeze(noisyExcitationResponseInstances(instanceIndex,:,targetConeID )), 'ko', 'LineWidth', 1.0);
    set(gca, 'YLim', [0 max(noisyExcitationResponseInstances(:))]);
    if (row == 8)
        xlabel('time (seconds)');
    end
    if (col == 1) && (row == 8)
        ylabel(sprintf('excitations per\nintegration time'));
    end
    title(sprintf('instance %d', instanceIndex));
end
