% Demonstrate the oiSequence class.
%
% Description
%    An oiSequence describes a dynamic retinal image, essentially a retinal
%    image video.  The oiSequence is not a general video, but it applies to
%    the case in which there is one basic stimulus that is either mixed
%    with a background or whose contrast is scaled over time.  Eye
%    movements are also included.  This simplification enables us to
%    compute many psychophysical stimuli efficiently; but it is not
%    completely general.
%
% Notes:
%    * For some purposes the oisCreate() function presents a simpler
%      approach to creating a sequence. That function has several built-in
%      cases, and more may be added.
%
% See Also: oisCreate
%

% History:
%    xx/xx/16  NPC  ISETBIO TEAM, 2016
%    10/19/18  JNM  Formatting

%% Generate a uniform scene
meanLuminance = 100;
uniformScene = sceneCreate('uniform equal photon', 128);
% square scene with desired FOV
FOV = 2.0;
uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
% 1 meter away
uniformScene = sceneSet(uniformScene, 'distance', 1.0);
uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);

% Generate a gabor scene
gaborParams = struct('freq', 2, 'contrast', 1.0, 'ph', 0, 'ang',  0, ...
    'row', 128, 'col', 128, 'GaborFlag', false);
gaborScene = sceneCreate('harmonic', gaborParams);
gaborScene = sceneSet(gaborScene, 'wAngular', FOV);
gaborScene = sceneSet(gaborScene, 'distance', 1.0);
gaborScene = sceneAdjustLuminance(gaborScene, meanLuminance);

%% generate the stimulus modulation function
stimulusSamplingInterval = 60 / 1000;
oiTimeAxis = -0.6:stimulusSamplingInterval:0.6;
stimulusRampTau = 0.165;

% monophasic modulation function
modulationFunction1 = ...
    0.7 * exp(-0.5 * (oiTimeAxis / stimulusRampTau) .^ 2);
modulationFunction2 = -modulationFunction1;

% biphasic modulation function
modulationFunction3 = 0.7 ...
    * (exp(-0.5 * ((oiTimeAxis - 0.2) / stimulusRampTau) .^ 2) ...
    - exp(-0.5 * ((oiTimeAxis + 0.2) / stimulusRampTau) .^ 2));

% Default human optics
oi = oiCreate('human');
oi = oiCompute(oi, uniformScene);

% Compute the background and the modulated optical images
oiBackground = oiCompute(oi, uniformScene);
oiModulated = oiBackground;
oiModulatedGabor = oiCompute(oi, gaborScene);
modulationRegion.radiusInMicrons = 250;

%% Adding a background
% oiSequence object for computing a sequence of ois where the oiModulated
% (uniform field) is ADDED to the oiBackground over an 250 micron radius
% region using a monophasic modulation function
theOIsequence(1) = oiSequence(oiBackground, oiModulated, ...
    oiTimeAxis, modulationFunction1, 'composition', 'add', ...
    'modulationRegion', modulationRegion);
theOIsequence(1).visualize('movie illuminance');

%% Blending with a background
% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground using a monophasic
% modulation function
theOIsequence(2) = oiSequence(oiBackground, oiModulated, ...
    oiTimeAxis, modulationFunction2, 'composition', 'add', ...
    'modulationRegion', modulationRegion);
theOIsequence(2).visualize('movie illuminance');

%% Adding a grating
% This is a oiSequence object for computing a sequence of ois where the
% oiModulated (a grating) is ADDED with the oiBackground using a monophasic
% modulation function
theOIsequence(3) = oiSequence(oiBackground, oiModulatedGabor, ...
    oiTimeAxis, modulationFunction1, 'composition', 'add');
theOIsequence(3).visualize('movie illuminance');

%% Blending a grating
% This is a oiSequence object for computing a sequence of ois where the
% oiModulated (a grating) is BLENDED with the oiBackground using a
% monophasic modulation function
theOIsequence(4) = oiSequence(oiBackground, oiModulatedGabor, ...
    oiTimeAxis, modulationFunction1, 'composition', 'blend');
theOIsequence(4).visualize('movie illuminance');

%%
% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground  over an 250 micron radius
% using a biphasic modulation function
theOIsequence(5) = oiSequence(oiBackground, oiModulatedGabor, ...
    oiTimeAxis, modulationFunction3, 'composition', 'blend', ...
    'modulationRegion', modulationRegion);
theOIsequence(5).visualize('movie illuminance');

%%
% Render a montage of the 5th oisequence
theOIsequence(5).visualize('montage');
