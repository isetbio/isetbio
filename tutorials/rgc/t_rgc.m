% funcftion [scene, sceneRGB, oi, sensor] = movieCreate(varargin)


% (HJ) ISETBIO TEAM, 2014
% (JRG) modified 10/2015

%% Init
ieInit;

%% Compute a Gabor patch scene

% Set up Gabor stimulus using sceneCreate('harmonic',params)
fov = 0.6;
params.freq = 6; params.contrast = 1;
params.ph  = 0;  params.ang = 0;
params.row = 64; params.col = 64;
params.GaborFlag = 0.2; % standard deviation of the Gaussian window

% Set up scene, oi and sensor
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

% These parameters are for other stuff.
params.expTime = 0.01;
params.timeInterval = 0.01;
params.nSteps = 5;     % Number of stimulus frames

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

sensor = sensorSet(sensor, 'exp time', params.expTime); 
sensor = sensorSet(sensor, 'time interval', params.timeInterval); 

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
fprintf('Computing cone isomerization:    ');

% ieSessionSet('wait bar',true);
wFlag = ieSessionGet('wait bar');
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
for t = 1 : params.nSteps
    if wFlag, waitbar(t/params.nSteps,wbar); end
        
    % Update the phase of the Gabor
    params.ph = 2*pi*(t-1)/params.nSteps; % one period over nSteps
    scene = sceneCreate('harmonic', params);
    scene = sceneSet(scene, 'h fov', fov);
    
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

% Show the movie of volts
vcNewGraphWin;axis image
for ii=1:params.nSteps
    imagesc(volts(:,:,ii)); pause(1); 
end

%% Outer segment calculation

% % % Input = cone current
os = osCreate('biophys');
% 
% % Compute the photocurrent
os = osCompute(os, sensor);
% 
% % Plot the photocurrent for a pixel
% osPlot(os,sensor);

% Input = RGB
% os = osCreate('identity');
% os = osSet(os, 'rgbData', sceneRGB);

%% Build rgc

% rgc1 = rgcCreate('GLM', scene, sensor, os, 'right', 3.0, 180);

clear params
params.scene = scene; 
params.sensor = sensor; 
params.outersegment = os;
params.eyeSide = 'left'; 
params.eyeRadius = 9; 
params.eyeAngle = 90;
rgc1 = rgcCreate('GLM', params);

rgc1 = rgcCompute(rgc1, os);

% rgcPlot(rgc1, 'mosaic');
% rgcPlot(rgc1, 'linearResponse');
rgcPlot(rgc1, 'spikeResponse');
%% Build rgc response movie
%  https://youtu.be/R4YQCTZi7s8

% % osLinear
% rgcMovie(rgc1, sensor);

% % osIdentity
% rgcMovie(rgc1, os);