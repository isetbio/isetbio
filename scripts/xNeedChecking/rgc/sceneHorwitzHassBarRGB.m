function [sceneRGB] = sceneHorwitzHassBarRGB(params)
% Build a scene object that consists of a movie of a block flashing on and
% off over a gray background. This is used for probing the RGC array.
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
% % sceneRGB = sceneHorwitzHassBarRGB(params);
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

barwidth = 5;

barMovie = zeros(params.image_size,params.image_size,3);
sceneRGB = zeros(params.image_size,params.image_size,params.nsteps,3);
fprintf('Computing RGB scene data    ');
for fstep = 1:params.nsteps

    fprintf('\b\b\b%02d%%', round(100*fstep/params.nsteps));

    barMovie = 0.5*ones(params.image_size,params.image_size,3);
%     barMovie(26:30,1+3*fstep+1:1+3*fstep+barwidth,:) = 0.5 + 0.499*ones(5,barwidth,3);% - mod(fstep,8)*0.999*ones(15,barwidth,3);
    barMovie(:,1+3*fstep+1:1+3*fstep+barwidth,:) = 0.5 + 0.499*ones(params.image_size,barwidth,3);% - mod(fstep,8)*0.999*ones(15,barwidth,3);
%     barMovie(:,17+1*fstep+1:17+1*fstep+barwidth,:) = 0.5 + 0.499*ones(params.image_size,barwidth,3);%
%     barMovie(17+1*fstep+1:17+1*fstep+barwidth,:,:) = 0.5 + 0.499*ones(barwidth,params.image_size,3);%
%     barMovie(:,17+1*fstep+1:17+1*fstep+barwidth,:) = 0.5 + 0.499*ones(params.image_size,barwidth,3) - mod(fstep,2)*0.999*ones(params.image_size,barwidth,3);
    
    barMovie(1,1,:) = 1;

    % % % % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(barMovie, 'rgb', meanLuminance, display);
    
    % sceneRGB(:,:,fstep,:) = sceneGet(scene,'rgb');
    sceneRGB(:,:,fstep,:) = barMovie;
    sceneRGB(1,1,fstep,:) = sceneRGB(2,2,fstep,:);
end
fprintf('     \n');

