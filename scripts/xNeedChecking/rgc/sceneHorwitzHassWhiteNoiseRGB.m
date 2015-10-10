function [sceneRGB] = sceneHorwitzHassWhiteNoiseRGB(params)
% Build a scene object that consists of an RGB white noise movie. 
% This is used for probing the RGC array.
%  
% Inputs: params structure consisting of image_size, meanLuminance,
% nsteps and fov.
% 
% Outputs: sceneRGB is a (image_size, image_size, nframes, 3) array, where
% the last dimension is RGB index.
% 
% Example:
% 
% params.image_size = 64;
% params.meanLuminance = 100;
% params.nsteps = 30;
% params.fov = 0.8;
% [scene, display] = sceneHorwitzHassWhiteNoise(params);
% oi  = oiCreate('wvf human');
% sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
% sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);
% % % % % % % % % % % % % % 


% Initialize isetBIO
s_initISET;

% show_movie_flag = params.disp_movie;

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
meanLuminance = params.meanLuminance;

% stimulusRGBdata = rgbGaborColorOpponentNormalized(params);
stimRGBraw = 0.5+(0.25*randn(params.image_size,params.image_size,3));
stimulusRGBdata = floor(254*abs(stimRGBraw)./max(stimRGBraw(:)));

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


fprintf('Computing RGB scene data    ');
for step = 1:params.nsteps
%     fprintf('     \n');
    fprintf('\b\b\b%02d%%', round(100*step/params.nsteps));
%     fprintf('     \n');
    
    stimRGBraw = 0.5+(0.25*randn(params.image_size,params.image_size,3));
    stimulusRGBdata = floor(254*abs(stimRGBraw)./max(stimRGBraw(:)));
    % % % % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(stimulusRGBdata, 'rgb', meanLuminance, display);
    
    sceneRGB(:,:,step,:) = sceneGet(scene,'rgb');
    
end
fprintf('     \n');


