% t_outersegment
%
% James Golden, 6/24/15
%
% Model cone responses for experiment on detection of gabor color opponents
% found in "Chromatic detection from cone photoreceptors to V1 neurons to
% behavior in rhesus monkeys" by Horwitz, Hass, Angueyra, Lindbloom-Brown &
% Rieke, J. Neuronscience, 2015
% 
% 

% clear classes
clear

% parameters found in Fig. 6 caption
params.color_val = 1;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
params.contrast = 0.8;        % set max contrast of gabor
params.image_size = 64;       % scene is image_size x image_size pixels
params.disp_movie = 1;        % display movie flag


params.fov = 1.2;                % degrees, sd = 0.4 deg, truncated at 3 sd
params.freq = 3.6;       % cyc/image = cyc/1.2 degs; want 3 cyc/deg = x cyc/1.2 degs, x = 3.6
t = 1;
params.period = 10; % in millseconds
params.ph  = 2*pi*((t-1)/params.period);  % vary with time at 3 hz?
params.ang = 0;
params.row = params.image_size; 
params.col = params.image_size;
params.GaborFlag = (128/params.image_size)*0.1/3.6; % standard deviation of the Gaussian window FIX
params.nsteps = 8;% total length of movie in millseconds
params.meanLuminance = 100;

% build scene
[scene, display] = build_scene_horwitz_hass_2015(params);

% build optical image
oi  = oiCreate('wvf human');

% build sensor
% Gabor color-opponent present
sensor = build_sensor_horwitz_hass_2015(params, scene, oi, display);

% % load sensor
% load('sensor_cv1_contp25.mat','sensor');

%%
% build outersegment
noise_flag = 1;
os1 = outersegment(noise_flag);

% compute outersegment output
tic
os1 = outersegmentCompute(os1,sensor,'linearfilter')
% osnl1 = outersegmentCompute(os1,'nonlinear')
toc

% plot signal and signal+noise from one cone 
[sz1 sz2 sz3] = size(os1.ConeCurrentSignal);
sig(1,:)          = os1.ConeCurrentSignal(round(sz1/2),round(sz2/2),:);
sigplusnoise(1,:) = os1.ConeCurrentSignalPlusNoise(round(sz1/2),round(sz2/2),:);

% figure; 
% plot(sig);
% hold on;
% plot(sigplusnoise,'r');


ConeCurrentSignalCell = convert_array_cell(params, sensor, os1);

% plot representation in LMS contrast space with max at each point in time
cs1 = max(ConeCurrentSignalCell{1});
cs2 = max(ConeCurrentSignalCell{2});
cs3 = max(ConeCurrentSignalCell{3});

% figure; scatter3(cs1(:),cs2(:),cs3(:),'b');

%% Find max values, compute 100 noisey versions 


reps_max = 100;
coneSamplingRate = 825;
cone_mosaic = sensorGet(sensor,'cone type');
[noisyL noisyM noisyS] = noisy_repeats(os1.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);

% % Running the function takes a couple minutes, so load instead
% load('noisyLMS.mat','noisyL','noisyM','noisyS');

%% Gabor color-opponent absent
clear sensor

% set contrast to 0%
params.contrast = 0;
sensor = build_sensor_horwitz_hass_2015(params, scene, oi, display);

% load('sensor_cv1_contp0.mat','sensor');

% build outersegment
noise_flag = 1;
os2 = outersegment(noise_flag,0);

% compute outersegment output
tic
os2 = outersegmentCompute(os2,sensor,'linearfilter')
% osnl1 = outersegmentCompute(os1,'nonlinear')
toc


ConeCurrentSignalCell2 = convert_array_cell(params, sensor, os1);

% plot representation in LMS contrast space with max at each point in time
csg1 = max(ConeCurrentSignalCell2{1});
csg2 = max(ConeCurrentSignalCell2{2});
csg3 = max(ConeCurrentSignalCell2{3});

% figure; scatter3(csg1(:),csg2(:),csg3(:),'b');

%% Compute 100 noisey version of cont = 0%

reps_max = 100;
coneSamplingRate = 825;
cone_mosaic = sensorGet(sensor,'cone type');
[noisyL_cont0 noisyM_cont0 noisyS_cont0] = noisy_repeats(os2.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);

% Running the function takes a couple minutes, so load instead
% load('noisyLMS_cont0.mat','noisyL_cont0','noisyM_cont0','noisyS_cont0');

%% LDA to find threshold

% half of input data is used as training set, half is used as test set
compute_lda(noisyL, noisyM, noisyS, noisyL_cont0, noisyM_cont0, noisyS_cont0);

