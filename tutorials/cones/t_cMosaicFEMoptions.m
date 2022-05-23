% Advanced usage of the new @cMosaic object with eye movements
%
% Description:
%    Shows how to change fEM properties with the new cone mosaic class,
%    @cMosaic, and different ways to center the generated fEM paths.
%
% See Also:
%   t_cMosaicSingleEyeMovementPaths
%   t_cMosaicMultipleEyeMovementPaths

% History:
%    09/20/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


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

%% Generate 1 eye movement path lasting for 400 msec
eyeMovementDurationSeconds = 400/1000;

% Change some fEM parameter (increase the positional noise)
defaultFEM = fixationalEM();
driftModelPositionNoiseStd = 1.2 * defaultFEM.positionNoiseStd;


t1Msec = 300;
t2Msec = 320;  

nTrials = 10;
randomSeed = 12354;

hFig = figure(1);
ax = subplot(1,2,1);
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'centerPaths', true, ...
        'driftModelPositionNoiseStd', driftModelPositionNoiseStd, ...
        'nTrials', nTrials, ...
        'randomSeed', randomSeed);

temporalSupportMsec = 1000*cm.fixEMobj.timeAxis;
[~,timePoint1] = min(abs(temporalSupportMsec-t1Msec));
[~,timePoint2] = min(abs(temporalSupportMsec-t2Msec));
t1Msec = temporalSupportMsec(timePoint1);
t2Msec = temporalSupportMsec(timePoint2);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'crossHairsOnMosaicCenter', true, ...
    'displayedEyeMovementData', struct('trial', 1:nTrials, 'timePoints', timePoint1:timePoint2), ...
    'plotTitle', sprintf('fEM data: %2.0f - %2.0f msec\n(centered paths)', t1Msec, t2Msec));

ax = subplot(1,2,2);
% Center fEM paths at t = 300 msec
centerFEMPathAtMsec = 300;
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'driftModelPositionNoiseStd', driftModelPositionNoiseStd, ...
        'centerPathsAtSpecificTimeMsec', centerFEMPathAtMsec, ...
        'nTrials', nTrials, ...
        'randomSeed', randomSeed);

cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'crossHairsOnMosaicCenter', true, ...
    'displayedEyeMovementData', struct('trial', 1:nTrials, 'timePoints', timePoint1:timePoint2), ...
    'plotTitle', sprintf('fEM data: %2.0f - %2.0f msec\n(paths centered at %2.0f msec)', t1Msec, t2Msec, centerFEMPathAtMsec));

