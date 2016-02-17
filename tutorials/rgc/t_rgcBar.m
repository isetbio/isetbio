%% t_rgcBar
% 
% A tutorial for the isetbio RGC object. A moving bar is created 
% in isetbio as a scene, and the sensor, outer segment and RGC objects are
% computed in response to the scene. 
% 
% (JRG) modified 10/2015

%% Init
ieInit;

%% Compute a Gabor patch scene

% Set up Gabor stimulus using sceneCreate('harmonic',params)
fov = 0.6;

params.barwidth = 5;
params.meanLuminance = 200;
params.row = 64; params.col = 64;
% params.freq = 6; params.contrast = 1;
% % params.ph  = 0;  params.ang = 0;
% params.row = 64; params.col = 64;
% params.GaborFlag = 0.2; % standard deviation of the Gaussian window

% Create display
display = displayCreate('CRT-Sony-HorwitzLab');

% Set up scene, oi and sensor
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

% These parameters are for other stuff.
params.expTime = 0.01;
params.timeInterval = 0.01;
params.nSteps = 50;     % Number of stimulus frames

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
fprintf('Computing cone isomerization:    \n');

% ieSessionSet('wait bar',true);
wFlag = ieSessionGet('wait bar');
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
for t = 1 : params.nSteps
    if wFlag, waitbar(t/params.nSteps,wbar); end
        

%     if t == 1 
        barMovie = zeros(params.row,params.col,3);
%     end
    barMovie(:,1+t:t+params.barwidth,:) = 0.5 + 0.499*ones(params.row,params.barwidth,3);
    barMovie(1,1,:) = 0;

    % % % % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);

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

%% Movie of the cone absorptions over cone mosaic

% coneImageActivity(sensor,'step',1,'dFlag',true);
%% Outer segment calculation

% Input = RGB
os = osCreate('identity');

coneSpacing = sensorGet(sensor,'width','um');
coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'coneSpacing', coneSpacing);

coneSampling = sensorGet(sensor,'time interval','sec');
os = osSet(os, 'coneSampling', coneSampling);

os = osSet(os, 'rgbData', sceneRGB);
% os = osCompute(sensor);

% % Plot the photocurrent for a pixel
% osPlot(os,sensor);
%% Build rgc

clear params

params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 12;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina = irCreate(os, params);

innerRetina.mosaicCreate('model','glm','mosaicType','on parasol');
%% Compute RGC response

innerRetina = irCompute(innerRetina, os);
for numberTrials = 1:10
    innerRetina.spikeCompute(innerRetina, os);
end

%%
% irPlot(innerRetina, 'mosaic');
% irPlot(innerRetina, 'linearResponse');
irPlot(innerRetina, 'rasterResponse');
% irPlot(innerRetina, 'psthResponse');

%% Build rgc response movie
% irMovie(rgc1, os);
%  https://youtu.be/R4YQCTZi7s8
