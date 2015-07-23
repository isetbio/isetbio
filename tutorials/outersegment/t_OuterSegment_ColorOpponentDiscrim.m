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

    params = paramsGaborColorOpponent();
    params.color_val = 1;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
    params.disp_movie = 0;        % display movie flag
        
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
    linearOS = osLinear('noiseflag',noiseflag);
    
    % compute outersegment current output
    linearOS.osCompute(sensor,params);
        
    % Find max values, compute 100 noisy versions    
    reps_max = 100;
    coneSamplingRate = 825;
    cone_mosaic = sensorGet(sensor,'cone type');
    [noisyL noisyM noisyS] = coneSignalNoisyRepeats(linearOS.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
    
    % % Running the function takes a couple minutes, so load instead
    % load('noisyLMS.mat','noisyL','noisyM','noisyS');
    
%% instantiate a LinearOuterSegment class for a zero-contrast gabor
    % build scene 
    clear scene display oi sensor
    params.contrast = 0;
    
    [scene, display] = sceneHorwitzHass(params);
    displayClose;
    
    % build optical image
    oi  = oiCreate('wvf human');
    
    % build sensor
    sensor = sensorHorwitzHass(params, scene, oi, display);
        
    % build outersegment
    noiseflag = 1;
    linearOS2 = osLinear('noiseflag',noiseflag);
    
    % compute outersegment current output
    linearOS2.osCompute(sensor,params);
        
    % Find max values, compute 100 noisy versions    
    reps_max = 100;
    coneSamplingRate = 825;
    cone_mosaic = sensorGet(sensor,'cone type');
    [noisyL_cont0 noisyM_cont0 noisyS_cont0] = coneSignalNoisyRepeats(linearOS2.ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic);
    
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