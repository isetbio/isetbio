%% cMosaic responses to a dead-leaves scene
%
% Compute the responses of a small cMosaic to a dead-leaves scene during a
% short fixational eye-movement sequence. The response at the final time
% point and the eye-movement path are shown in the cMosaic window.
%
% See also
%   cMosaic, sceneDeadleaves, cMosaic.compute, cMosaic.emGenSequence

%%
ieInit;

%% Create the dead-leaves scene and optical image

rng(1);
sceneSize = 128;
sceneSigma = 3;
sceneFOV = 0.6;

scene = sceneDeadleaves(sceneSize, sceneSigma);
scene = sceneSet(scene, 'fov', sceneFOV);

oi = oiCreate('human');
oi = oiCompute(oi, scene, 'pad value', 'mean');

%% Load a small cMosaic and create a short eye-movement sequence

cm = mosaicLoad([0.5 0.5], [0 0]);
cm.integrationTime = 10/1000;
cm.noiseFlag = 'none';

eyeMovementDurationSeconds = 100/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
    'microsaccadeType', 'none', ...
    'nTrials', 1, ...
    'randomSeed', 1);

%% Compute the noise-free cone excitations

[excitations, ~, ~, ~, timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true);

%% Show the final response and the eye-movement path

cm.plot('excitations and eye movements', excitations, ...
    'time point', 'last', ...
    'plot title', sprintf('Dead leaves response at %d ms', ...
        round(timeAxis(end)*1000)));

%%
