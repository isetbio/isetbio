% t_fitConeOutputWithGLM
% 
% The current computational model of RGC response maps the raw pixel
% stimuli directly to RGC spikes. We want to do a new type of fit to the
% RGC spikes using simulated cone currents as the input stimulus, as
% opposed to the raw pixel data. In order to do this, we must simulate the
% cone current responses to the movie used to fit the RGC spikes in the
% lab. 
% 
% Here, we initialize isetbio scene, optical image, sensor (cone
% isomerizations) and cone outer segment objects and compute the cone outer
% segment current in response to the binary white noise movie used to fit
% RGC responses. We use a particular movie from an experiment on piece
% 2013-08-19-6 in the Chichilnisky lab and load the RGC spikes from 

%% Init
clear
ieInit;

%% Compute a Gabor patch scene

% Set up Gabor stimulus using sceneCreate('harmonic',params)
params.fov = 3; % degrees

params.meanLuminance = 200;
params.row = 80; params.col = 40;
% params.freq = 6; params.contrast = 1;
% % params.ph  = 0;  params.ang = 0;
% params.row = 64; params.col = 64;
% params.GaborFlag = 0.2; % standard deviation of the Gaussian window

% Create display
display = displayCreate('CRT-Sony-HorwitzLab');

% Set up scene, oi and sensor
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', params.fov);
% vcAddObject(scene); sceneWindow;

% These parameters are for other stuff.
params.expTime      = 8.3275e-03/2;
params.timeInterval = 8.3275e-03/2;
params.nSteps = 20;     % Number of stimulus frames

%% Initialize the optics and the sensor


oi  = oiCreate('wvf human');

coneP = coneCreate; % The cone properties properties

% Check out experiment folder to get these exact numbers
retinalRadiusDegrees = 40;
retinalPolarDegrees  = 0;
retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; whichEye = 'right';
sensor = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
sensor = sensorSetSizeToFOV(sensor, params.fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', params.expTime); % 1 ms
sensor = sensorSet(sensor, 'time interval', params.timeInterval); % 1 ms

% This computes with no sensor or photon noise, just the mean
sensor = sensorSet(sensor,'noise flag',0);


% % % % % % 
% sensor = sensorCreate('human');
% sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
% 
% sensor = sensorSet(sensor, 'exp time', params.expTime); 
% sensor = sensorSet(sensor, 'time interval', params.timeInterval); 

%% Compute a dynamic set of cone absorptions

% We want to produce a scene video that translates into an oi video that
% becomes a cone absorption video.  At present coneAbsorptions ONLY does
% this using eye movements, not by creating a series of images.  This code
% represents our first effort to produce dynamic scenes.

% We are literally going to recreate a set of scenes with different phase
% positions and produce the scenes, ois, and cone absorptions by the loop.
% The result will be a time series of the cone photon absorptions.
%
% We are reluctant to make scene(:,:,:,t) because we are frightened about
% the size.  But it still might be the right thing to do.  So the code here
% is an experiment and we aren't sure how it will go.

% sceneRGB = zeros([sceneGet(scene, 'size') params.nSteps 3]); % 3 is for R, G, B
% sensorPhotons = zeros([sensorGet(sensor, 'size') params.nSteps]);
% stimulus = zeros(1, params.nSteps);
fprintf('Computing cone isomerization:    \n');

% ieSessionSet('wait bar',true);
wFlag = ieSessionGet('wait bar');
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Load binary white noise fit movie from Chichilnisky Lab
% We have spike outputs from RGCs in response to this movie
load('/Users/james/Documents/MATLAB/akheitman/fitmovie_8pix_Identity_8pix.mat')
% figure; imagesc(BWmovie.fitmovie.movie_byblock{1}.matrix(:,:,1))

% Loop through frames to build movie
for t = 1 : params.nSteps
    if wFlag, waitbar(t/params.nSteps,wbar); end
        

    if t == 1 
        stimulusRGBdata = zeros(params.row,params.col,1);
    end
   
    stimulusRGBdata = BWmovie.fitmovie.movie_byblock{1}.matrix(:,:,floor(t/2)+1);
    
    % % % % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);

    scene = sceneSet(scene, 'h fov', params.fov);

    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions
    sensor = sensorCompute(sensor, oi);

    if t == 1
        volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
    end
    
    volts(:,:,t) = sensorGet(sensor, 'volts');
    
    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
sensor = sensorSet(sensor, 'volts', volts);
% vcAddObject(sensor); sensorWindow;


%% Movie of the cone absorptions over cone mosaic

% coneImageActivity(sensor,'step',1,'dFlag',true);
%% Outer segment calculation

% Input = RGB
os = osCreate('biophys');

coneSpacing = sensorGet(sensor,'width','m');
os = osSet(os, 'patchsize', coneSpacing);

coneSampling = sensorGet(sensor,'time interval','sec');
os = osSet(os, 'timestep', coneSampling);

os = osCompute(os,sensor);
% % Plot the photocurrent for a pixel
osPlot(os,sensor);
