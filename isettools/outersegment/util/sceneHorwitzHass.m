function [scene , display] = sceneHorwitzHass(params)
%Build a scene object following parameters from "Chromatic detection from 
% cone photoreceptors to V1 neurons to behavior in rhesus monkeys" by 
% Horwitz, Hass, Angueyra, Lindbloom-Brown & Rieke, J. Neuronscience, 2015
%  
% % build scene
% % parameters found in Fig. 6 caption
% params.color_val = 1;         % 1-4, choose color opponent gabor
% params.contrast_set = 0.8;    % set max contrast of gabor
% params.image_size = 64;       % scene is image_size x image_size pixels
% params.disp_movie = 0;        % display movie flag
% % 
% fov = 1.2; % degrees, sd = 0.4 deg, truncated at 3 sd
% params.color_val = cv; % 1 = s_iso
% params.freq = 3.6; % cyc/image = cyc/1.2 degs; want 3 cyc/deg = x cyc/1.2 degs, x = 3.6
% params.contrast = cont_mult;
% t = 0; params.ph  = 2*pi*((t-1)/333);  % vary with time at 3 hz?
% params.ang = 0;
% params.row = 128; params.col = 128;
% params.GaborFlag = 0.1/3.6; % standard deviation of the Gaussian window FIX
% 
% scene = build_scene_horwitz_hass_2015(params);
% % % % % % % % % % % % % % 


% Initialize isetBIO
s_initISET;

show_movie_flag = params.disp_movie;

% Display information you need to get from Rieke/Horwitz
NspectralSamples = 81;
%%%% For some reason background gray level is changing
load('spdCalibrationHorwitzHass.mat');
% wave  = zeros(NspectralSamples, 1);  % e.g. 400 420 440 .. 800
wave = cal.monSpectWavelengths; 
spd = reshape(cal.monSpect,81,3); % power at each of the above wavelengths, in Watts/steradian/m^2/nm, for each of the R,G,B channels
dpi   = 96;              % display resolution in pixels per inch
vd_inMeters = 2.0;       % viewing distance in meters

% Stimulus information

cv = params.color_val;
cont_mult = params.contrast;

% clear sensor oi volts
clear volts

% % parameters found in Fig. 6 caption of Hass, et al., 2015
% fov = 1.2; % degrees, sd = 0.4 deg, truncated at 3 sd
% params.color_val = cv; % 1 = s_iso
% params.freq = 3.6; % cyc/image = cyc/1.2 degs; want 3 cyc/deg = x cyc/1.2 degs, x = 3.6
% params.contrast = cont_mult;
% t = 0; params.ph  = 2*pi*((t-1)/333);  % vary with time at 3 hz?
% params.ang = 0;
% % params.row = image_size; params.col = image_size;
% params.GaborFlag = (128/params.image_size)*0.1/3.6; % standard deviation of the Gaussian window FIX

% Stimulus information
meanLuminance = params.meanLuminance;

stimulusRGBdata = rgbGaborColorOpponentNormalized(params);
% figure; hist(sqrt(sum(reshape(stimulusRGBdata.^2,params.row*params.row,3),2)),40)

if show_movie_flag == 1
    playGaborMovie(params.color_val,params);
    close;
end

% load('monitor_cal.mat');

% Generate a display object to model Horwitz's display
display = displayCreate;


% Set the display's name
display = displaySet(display, 'name', 'Horwitz');


% Set the display's SPDs
display = displaySet(display, 'wave', wave);
display = displaySet(display, 'spd', spd);
display = displaySet(display, 'ambientspd', zeros(1,length(spd))); 

% Set the display's resolution (dots-per-inch)
display = displaySet(display, 'dpi', dpi);


% Set the display's viewing distance
display = displaySet(display, 'viewing distance', vd_inMeters);


% % % % Generate scene object from stimulus RGB matrix and display object
scene = sceneFromFile(stimulusRGBdata, 'rgb', meanLuminance, display);
vcAddAndSelectObject(scene);
sceneWindow

