% James Golden & Nicolas Cottaris, 7/14/15
%
% Model cone responses for experiment on detection of gabor color opponents
% found in "Chromatic detection from cone photoreceptors to V1 neurons to
% behavior in rhesus monkeys" by Horwitz, Hass, Angueyra, Lindbloom-Brown &
% Rieke, J. Neuroscience, 2015
%
%
%%  preliminaries - define gabor color opponent stimulus, scene, oi, sensor

clear

% parameters found in Fig. 6 caption

params = paramsGaborColorOpponent()
params.color_val = 2;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
params.disp_movie = 1;        % display movie flag

% build scene
[scene, display] = sceneHorwitzHass(params);
displayClose;
% build optical image
oi  = oiCreate('wvf human');

% build sensor
% Gabor color-opponent present
sensor = sensorHorwitzHass(params, scene, oi, display);

% % load sensor
% load('sensor_cv1_contp25.mat','sensor');

%% instantiate a LinearOuterSegment class for a high-contrast gabor

% build outersegment
noiseflag = 1;
linearOS = osLinear('noiseFlag',noiseflag);
% adaptedOS = osBioPhys('noiseFlag',noiseflag);
% compute outersegment current output
linearOS = osLinearCompute(linearOS,sensor,params);
% adaptedOS = osBioPhysCompute(adaptedOS,sensor,params);
% Find max values, compute 100 noisy versions
reps_max = 100;
coneSamplingRate = 825;
cone_mosaic = sensorGet(sensor,'cone type');
[noisyL noisyM noisyS] = coneSignalNoisyRepeats(linearOS.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);

% [sz1, sz2, sz3] = size(adaptedOS.ConeCurrentSignal);
% adaptedOS_ConeCurrent = adaptedOS.ConeCurrentSignal - repmat(adaptedOS.ConeCurrentSignal(:,:,sz3),[1,1,sz3]);
% % [noisyL noisyM noisyS] = coneSignalNoisyRepeats_allpts(adaptedOS_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);
% % [noisyL noisyM noisyS] = coneSignalNoisyRepeats(adaptedOS_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);
% % [noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats_allpts(linearOS2.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
% % [noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats_allpts(adaptedOS_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);
% [noisyContrastFull] = coneSignalNoisyRepeats_sp(adaptedOS_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);

% % Running the function takes a couple minutes, so load instead
% load('noisyLMS.mat','noisyL','noisyM','noisyS');

%% instantiate a LinearOuterSegment class for a zero-contrast gabor
% build scene
% clear scene display oi sensor
params.contrast = 0;

[scene2, display2] = sceneHorwitzHass(params);
displayClose;

% build optical image
oi2  = oiCreate('wvf human');

% build sensor
sensor2 = sensorHorwitzHass(params, scene2, oi2, display2);

% build outersegment
noiseflag = 1;
linearOS2 = osLinear('noiseFlag',noiseflag);

% adaptedOS2 = osBioPhys('noiseFlag',noiseflag);
% compute outersegment current output
linearOS2 = osLinearCompute(linearOS2,sensor2,params);

% adaptedOS2 = osBioPhysCompute(adaptedOS2,sensor2,params);

% adaptedOS2_ConeCurrent = adaptedOS2.ConeCurrentSignal - repmat(adaptedOS2.ConeCurrentSignal(:,:,sz3),[1,1,sz3]);

% Find max values, compute 100 noisy versions
reps_max = 100;
coneSamplingRate = 825;
cone_mosaic = sensorGet(sensor2,'cone type');
[noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats(linearOS2.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
% [noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats(adaptedOS2_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);
% [noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats_allpts(adaptedOS2_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);
% [noisyContrastZero] = coneSignalNoisyRepeats_sp(adaptedOS2_ConeCurrent, reps_max, coneSamplingRate, cone_mosaic);

% % Running the function takes a couple minutes, so load instead
% load('noisyLMS_cont0.mat','noisyL_cont0','noisyM_cont0','noisyS_cont0')

%% LDA to find threshold

% half of data is used as training set, half is used as test set
classifyLinearDiscr(noisyL, noisyM, noisyS, noisyL_cont0, noisyM_cont0, noisyS_cont0);

%%
%     % instantiate a NonLinearOuterSegment class
%     noiseflag = 1;
%     nonlinearOS = NonLinearOuterSegment('noiseflag',noiseflag);
%     nonlinearOS.temporalFilter(sensor,params);
%