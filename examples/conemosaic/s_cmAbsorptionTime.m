%% Check cMosaic excitation dimensions across time
%
% cMosaic.compute returns excitations in trials-by-time-by-cones
% format. Unlike the legacy rectangular cone mosaic, a cMosaic is an
% irregular lattice, so its responses are not stored as
% rows-by-columns-by-time. This example checks single- and
% multi-time-point computations and shows how to select a
% one-dimensional sample of 100 cones.
%
% See also
%   cMosaic, cMosaic.compute, cMosaic.emGenSequence, cMosaic.plot
%
% BW, ISETBIO Team, 2016. Modernized 2026.

%% Initialize a uniform scene and human optical image

ieInit;

scene = sceneCreate('uniform ee', 256);
scene = sceneSet(scene, 'fov', 0.5);
oi = oiCreate('human');
oi = oiCompute(oi, scene, 'pad value', 'mean');

%% A two-dimensional mosaic with one temporal sample

cm = cMosaic( ...
    'sizeDegs', [0.25 0.25], ...
    'integrationTime', 0.05, ...
    'noiseFlag', 'none', ...
    'useParfor', false);
cm.emGenSequence(cm.integrationTime, 'microsaccadeType', 'none', ...
    'randomSeed', 1);
excitations = cm.compute(oi, 'withFixationalEyeMovements', true);

% excitations is [N Trials, N Eye Positions, N cones]
assert(isequal(size(excitations), [1 1 cm.conesNum]));
cm.plot('excitations', excitations, ...
    'plot title', 'One temporal sample');

%% The same mosaic with three temporal samples

cm.integrationTime = 0.01;
nTimePoints = 3;
cm.emGenSequence(nTimePoints*cm.integrationTime, ...
    'microsaccadeType', 'none', 'randomSeed', 1);
excitations = cm.compute(oi, 'withFixationalEyeMovements', true);

assert(isequal(size(excitations), [1 nTimePoints cm.conesNum]));
cm.plot('excitations', excitations, 'time point', 'last', ...
    'plot title', 'Final response of a three-sample sequence');

%% A one-dimensional sample of 100 cones across the mosaic

% The modern cMosaic does not have a 1-by-100 rectangular-array
% constructor. Instead, select 100 cones nearest the horizontal
% meridian and retain their responses as a vector indexed by cone.
nSampledCones = 100;

% coneRFpositionsDegs(:,2) is each cone's vertical (y) position.
% Sorting abs(y) puts cones closest to y = 0—the horizontal meridian
[~, coneOrder] = sort(abs(cm.coneRFpositionsDegs(:, 2)));
sampledConeIndices = coneOrder(1:nSampledCones);
oneDimensionalExcitations = excitations(:, :, sampledConeIndices);

assert(isequal(size(oneDimensionalExcitations), ...
    [1 nTimePoints nSampledCones]));

ieFigure;
plot(cm.coneRFpositionsDegs(sampledConeIndices, 1), ...
    squeeze(oneDimensionalExcitations(1, end, :)), 'k.');
xlabel('Horizontal position (deg)');
ylabel('Excitations');
title('Excitations of 100 cones nearest the horizontal meridian');
grid on;

%% A two-dimensional mosaic with 20 temporal samples

nTimePoints = 20;
cm.emGenSequence(nTimePoints*cm.integrationTime, ...
    'microsaccadeType', 'none', 'randomSeed', 1);
excitations = cm.compute(oi, 'withFixationalEyeMovements', true);

assert(isequal(size(excitations), [1 nTimePoints cm.conesNum]));
cm.plot('excitations', excitations, 'time point', 'last', ...
    'plot title', 'Response of a 20-sample sequence');

%%
