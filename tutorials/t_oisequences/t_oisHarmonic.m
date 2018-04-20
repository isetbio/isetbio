%% t_oisHarmonic
% Demonstrates how to use the oiSequence class to generate an oiSequence of
% a sinusoid whose contrast is ramped on and off using a Gaussian temporal 
% modulation function.
%
% See also: oisCreate
%
% ISETBIO TEAM, 2017

%% 
ieInit

%% Create the default oi
oi = oiCreate('wvf human');

%% Spatial parameters
imgFov = .5;  % image field of view
vDist = 0.3;  % viewing distance (meter)

%% The stimulus
params.freq = 10;       % spatial frequencies of 10
params.contrast = 0.6;  % 60% contrast
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'name', sprintf('F %d', params.freq));
scene = sceneSet(scene, 'h fov', imgFov);
scene = sceneSet(scene, 'distance', vDist);
oiModulated =  oiCompute(oi, scene);
ieAddObject(scene);


%% The background
clear params   
params.contrast = 0;  % zero contrast
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'name', 'Background');
scene = sceneSet(scene, 'h fov', imgFov);
scene = sceneSet(scene, 'distance', vDist);
oiBackground =  oiCompute(oi, scene);
ieAddObject(scene);

%% Build the oiSequence
% We build the stimulus using a time series of weights. We have the mean
% field on for a while, then rise/fall, then mean field.
zTime = 50;  % Mean field beginning and end (ms)
stimWeights = fspecial('gaussian', [1, zTime], 15);
stimWeights = ieScale(stimWeights, 0, 1);
weights = [zeros(1, zTime), stimWeights, zeros(1, zTime)];

%% Temporal samples.  Typically 1 ms, which is set by the parameter in the
% cone mosasic integration time.  That time is locked to the eye movements.
tSamples = length(weights);
sampleTimes = 0.002 * (1:tSamples);  % Time in sec

% vcNewGraphWin;  plot(1:tSamples, weights, 'o');
% xlabel('Time (ms)');

% The weights define some amount of the constant background and some amount
% of the line on the same constant background
oiHarmonicSeq = oiSequence(oiBackground, oiModulated, ...
    sampleTimes, weights, 'composition', 'blend');
oiHarmonicSeq.visualize('movie illuminance');
