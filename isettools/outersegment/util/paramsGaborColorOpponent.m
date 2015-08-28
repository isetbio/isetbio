function params = paramsGaborColorOpponent()
% Model cone responses for experiment on detection of gabor color opponents
% found in "Chromatic detection from cone photoreceptors to V1 neurons to
% behavior in rhesus monkeys" by Horwitz, Hass, Angueyra, Lindbloom-Brown &
% Rieke, J. Neuroscience, 2015
% 
% 
    % parameters found in Fig. 6 caption
    params.color_val = 2;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
    params.contrast = 0.01;        % set max contrast of gabor
    params.image_size = 64;       % scene is (image_size X image_size) pixels
    params.disp_movie = 0;        % display movie flag
    
    params.fov = 1.2;                % degrees, sd = 0.4 deg, truncated at 3 sd
    params.freq = 3.6;       % cyc/image = cyc/1.2 degs; want 3 cyc/deg = x cyc/1.2 degs, x = 3.6
    params.period = 10; % in millseconds
    params.ph  = 2*pi*((0)/params.period);  % vary with time at 3 hz?
    params.ang = 0;
    params.row = params.image_size;
    params.col = params.image_size;
    params.GaborFlag = (128/params.image_size)*0.1/3.6; % standard deviation of the Gaussian window FIX
    params.nsteps = 666;% total length of movie in millseconds
    params.meanLuminance = 100;