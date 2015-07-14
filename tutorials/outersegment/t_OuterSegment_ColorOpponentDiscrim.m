% James Golden & Nicolas Cottaris, 7/14/15
%
% Model cone responses for experiment on detection of gabor color opponents
% found in "Chromatic detection from cone photoreceptors to V1 neurons to
% behavior in rhesus monkeys" by Horwitz, Hass, Angueyra, Lindbloom-Brown &
% Rieke, J. Neuroscience, 2015
% 
% 
%%  preliminaries - define gabor color opponent stimulus, scene, oi, sensor  
    % parameters found in Fig. 6 caption

    params = gabor_color_opponent_params();
    params.color_val = 1;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
    params.disp_movie = 0;        % display movie flag
        
    % build scene
    [scene, display] = build_scene_horwitz_hass_2015(params);
    displayClose;
    % build optical image
    oi  = oiCreate('wvf human');
    
    % build sensor
    % Gabor color-opponent present
    sensor = build_sensor_horwitz_hass_2015(params, scene, oi, display);
    
    % % load sensor
    % load('sensor_cv1_contp25.mat','sensor');
    
%% instantiate a LinearOuterSegment class for a high-contrast gabor
    
    % build outersegment
    noiseflag = 1;
    linearOS = LinearOuterSegment('noiseflag',noiseflag);
    
    % compute outersegment current output
    linearOS.temporalFilter(sensor,params);
        
    % Find max values, compute 100 noisy versions    
    reps_max = 100;
    coneSamplingRate = 825;
    cone_mosaic = sensorGet(sensor,'cone type');
    [noisyL noisyM noisyS] = noisy_repeats(linearOS.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
    
    % % Running the function takes a couple minutes, so load instead
    % load('noisyLMS.mat','noisyL','noisyM','noisyS');
    
%% instantiate a LinearOuterSegment class for a zero-contrast gabor
    % build scene 
    clear scene display oi sensor
    params.contrast = 0;
    
    [scene, display] = build_scene_horwitz_hass_2015(params);
    displayClose;
    
    % build optical image
    oi  = oiCreate('wvf human');
    
    % build sensor
    sensor = build_sensor_horwitz_hass_2015(params, scene, oi, display);
        
    % build outersegment
    noiseflag = 1;
    linearOS2 = LinearOuterSegment('noiseflag',noiseflag);
    
    % compute outersegment current output
    linearOS2.temporalFilter(sensor,params);
        
    % Find max values, compute 100 noisy versions    
    reps_max = 100;
    coneSamplingRate = 825;
    cone_mosaic = sensorGet(sensor,'cone type');
    [noisyL_cont0 noisyM_cont0 noisyS_cont0] = noisy_repeats(linearOS2.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
    
    % % Running the function takes a couple minutes, so load instead
    % load('noisyLMS_cont0.mat','noisyL_cont0','noisyM_cont0','noisyS_cont0')

    %% LDA to find threshold

    % half of data is used as training set, half is used as test set
    compute_lda(noisyL, noisyM, noisyS, noisyL_cont0, noisyM_cont0, noisyS_cont0);

%%
%     % instantiate a NonLinearOuterSegment class
%     noiseflag = 1;
%     nonlinearOS = NonLinearOuterSegment('noiseflag',noiseflag);
%     nonlinearOS.temporalFilter(sensor,params);
%     