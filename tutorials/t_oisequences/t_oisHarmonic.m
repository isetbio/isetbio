% Generate a sinusoid oiSequence with contrast ramping by Gaussian tmf.
%
% Description:
%    Demonstrates how to use the oiSequence class to generate an oiSequence
%    of a sinusoid whose contrast is ramped on and off using a Gaussian
%    temporal modulation function.
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
params.row = 256;       % 256 x 256 pixels
params.col = 256;
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'name', sprintf('F %d', params.freq));
scene = sceneSet(scene, 'h fov', imgFov);
scene = sceneSet(scene, 'distance', vDist);
oiModulated =  oiCompute(oi, scene);
ieAddObject(scene);

%% The background
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

%% Temporal samples.
% These are typically 1 ms, which is set by the parameter in the cone
% mosasic integration time.  That time is locked to the eye movements.
tSamples = length(weights);
sampleTimes = 0.002 * (1:tSamples);  % Time in sec
% vcNewGraphWin;
% plot(1:tSamples, weights, 'o');
% xlabel('Time (ms)');

%% Generate oiSequence
oiHarmonicSeq = oiSequence(oiBackground, oiModulated, ...
    sampleTimes, weights, 'composition', 'blend');

%% Visualize as a movie
oiHarmonicSeq.visualize('movie illuminance');

%% Visualize all frames in a montage
oiHarmonicSeq.visualize('montage');
