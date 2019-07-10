% Script demonstrating how to use the oiSequences function.
%
% Description:
%    The oisCreate function produces some classic psychophysical stimuli as
%    oiSequences. This script shows how to use that function.
%
%    The script creates Gaussian envelope harmonic (Gabor) oiSequences for
%    both monochrome and color, Vernier stimuli, and a flash.
%
%    Implementing directly with the oiSequence() class is illustrated in
%    the script t_oiSequence.m
%
% Notes:
%

% History:
%    xx/xx/17  BW   ISETBIO TEAM, 2017
%    10/18/18  JNM  Formatting

%%
ieInit

%% Harmonic - monochrome
% This code uses sceneCreate('harmonic', ...) function to create harmonics.
% We create a pair where one has positive contrast and a second has 0
% contrast but equal mean. This produces a contrast modulating oi sequence
clear hparams

% The modulating harmonic parameters. The possibilities are defined in the
% sceneCreate('harmonic', ...) function, and the complete list of
% parameters is returned by harmonicP.
hparams(2) = harmonicP;
hparams(2).freq = 6;
hparams(2).GaborFlag = .2;

% The matched, zero contrast, harmonic parameters
hparams(1) = hparams(2);
hparams(1).contrast = 0;

% The general scene properties can also be set. Here we make the scene one
% deg of visual angle
sparams.fov = 1;

% And then we make a Gaussian temporal modulation that brings the stimulus
% on and off
stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);

% The two harmonics are 'blended', which means at each moment in time we
% have a weighted sum of the two where the weights sum to one.
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);

vname = fullfile(isetbioRootPath, 'local', 'oisVideo1.mp4');
ois.visualize('movieilluminance', 'vname', vname);

% Save as a gif for the wiki page.
% gifName = fullfile(isetbioRootPath, 'wiki', 'images', 'oisHarmonic.gif');
% ieGIF(uData.movie, 'gifName', gifName);

% Now, show the time series of weights
ois.visualize('weights');

%% Change the spatial frequency parameter
hparams(2).freq = 10;
hparams(1) = hparams(2);
hparams(1).contrast = 0;

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, ...
    'sceneParameters', sparams);
ois.visualize('movie illuminance');

%% A little phase shifting, causing apparent motion

hparams(2).freq = 2;
hparams(2).ph =  pi / 6;
hparams(1) = hparams(2);
hparams(1).ph = 0;

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);

ois.visualize('movie illuminance');

%% A color Gabor patch

% Set up the color SPDs for background and test modulation
dsp = displayCreate('LCD-Apple.mat');
wave = displayGet(dsp, 'wave');
backSPD = displayGet(dsp, 'spd primaries') * 0.5 * ones(3, 1);
backSPD = Energy2Quanta(wave, backSPD);
[~, modSPD] = humanConeIsolating(dsp);
modSPD = Energy2Quanta(wave, modSPD);

clear hparams
hparams(2) = harmonicP;
hparams(2).freq = 6;
hparams(2).GaborFlag = 0.2;
hparams(2).ang = pi / 6;

hparams(2).backSPD = backSPD;
hparams(2).modSPD = modSPD(:, 1);  % Could be a mixture of the modSPDs
hparams(2).wave = wave;
hparams(1) = hparams(2);
hparams(1).contrast = 0;

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);

vname = fullfile(isetbioRootPath, 'local', 'oisVideo2.mp4');
ois.visualize('movie rgb', 'vname', vname);

%% Vernier
% The line offset case is managed using the sceneCreate('vernier', ...).
% [Note: BW - The code for the vernier stimuli needs simplification,
% checking, and improvement.]

clear vparams;
% The vernier parameter defaultss
vparams(2) = vernierP;

% A black background
vparams(2).name = 'offset';
vparams(2).bgColor = 0;

% The bar is added to the black background
vparams(1) = vparams(2);
vparams(1).barWidth = 0;
vparams(1).bgColor = 0.5;
vparams(1).name = 'uniform';

% The scene
sparams.fov = 1;

% The temporal modulation
stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);

% You can return the scenes, if you like.
[vernier, scenes] = oisCreate('vernier', 'add', stimWeights, ...
    'testParameters', vparams, 'sceneParameters', sparams);
vernier.visualize('movie illuminance');

% To see the scenes, prior to creating the optical image
% ieAddObject(scenes{1});
% ieAddObject(scenes{2});
% sceneWindow;

% Change up the parameters
vparams(2).barColor = 0.5;

vernier = oisCreate('vernier', 'add', stimWeights, ...
    'testParameters', vparams, 'sceneParameters', sparams);
vernier.visualize('movie illuminance');

%% Impulse (temporal)
clear iparams

sparams.fov = 1;
sparams.luminance = 100;
stimWeights = zeros(1, 50);
stimWeights(20:23) = 1;

impulse = oisCreate('impulse', 'add', stimWeights, ...
    'sceneParameters', sparams);

impulse.visualize('movie illuminance');
