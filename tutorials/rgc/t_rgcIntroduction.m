
%% t_rgcIntroduction
%
% In which our heroes lay out the basic architecture of the rgc class and
% its specializations.
%
% JG/BW ISETBIO Team, Copyright 2015
% (HJ) ISETBIO TEAM, 2014
% (JRG) modified 10/2015
%
% Aspirational:
%    function [scene, sceneRGB, oi, sensor] = movieCreate(varargin)
%

ieInit

%% Aspirational:
% We probably want something like
%   params = SET UP
%   sensor = sensorMovie('gabor','params',params,'oi',oi,'sensor',sensor);

% For realistic testing we need dynamic scenes.  We are in the processing
% of developing methods for automating the creation of dynamic scenes.  But
% this isn't yet done.  We are thinking we will store the dynamic scenes as
% movies in the sensor mosaic (cone absorptions).
%
% At some point we will want to have t_dynamicScene, t_dynamicOI,
% t_dynamicAbsorptions.  These will be cases in which the scene itself
% changes, or might even be a movie.
%
% We also have dynamics introduced by eye movements only.  That case is
% handled somewhat differently.
%
% For a case like this, there is no particular reason to store the dynamic
% scene.  If we are going to do a lot of experiments with the oi and sensor
% fixed, changing only the scene, then we need to recompute from scratch
% each time.  So we would loop on each scene to produce the video of sensor
% catches.

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

%% Compute a scene

% Set up Gabor stimulus using sceneCreate('harmonic',params)
% 
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
params.expTime = 0.005;
params.timeInterval = 0.005;
params.nSteps = 10;%60;     % Number of stimulus frames
params.nCycles = 4;
%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');
absorptions = sensorCreate('human');
absorptions = sensorSetSizeToFOV(absorptions, fov, scene, oi);

absorptions = sensorSet(absorptions, 'exp time', params.expTime); 
absorptions = sensorSet(absorptions, 'time interval', params.timeInterval); 

%% Compute a dynamic set of cone absorptions

wFlag = ieSessionGet('wait bar');
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
for t = 1 : params.nSteps
    if wFlag, waitbar(t/params.nSteps,wbar); end
        
    % Update the phase of the Gabor
    %
    params.ph = (2*pi)*params.nCycles*(t-1)/params.nSteps; % one period over nSteps
    scene = sceneCreate('harmonic', params);
    scene = sceneSet(scene, 'h fov', fov);
    
    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions
    absorptions = sensorSet(absorptions,'noise flag',0);
    absorptions = sensorCompute(absorptions, oi);

    if t == 1
        isomerizations = zeros([sensorGet(absorptions, 'size') params.nSteps]);
    end
    
    isomerizations(:,:,t) = sensorGet(absorptions, 'photons');
    
    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
absorptions = sensorSet(absorptions, 'photons', isomerizations);
% vcAddObject(sensor); sensorWindow;

% %% Movie of the cone absorptions over cone mosaic
% % from t_VernierCones by HM
% 
% step = 1;   % Step is something about time?
% % Display gamma preference could be sent in here
% tmp = coneImageActivity(absorptions,[],step,false);
% 
% % Show the movie
% % Aspirational
% % showMovie(tmp)
% 
% vcNewGraphWin;
% tmp = tmp/max(tmp(:));
% for ii=1:size(tmp,4)
%     img = squeeze(tmp(:,:,:,ii));
%     imshow(img); truesize;
%     title('Cone absorptions')
%     drawnow
% end

%% Outer segment calculation

% The outer segment converts cone absorptions into cone photocurrent.
% There are 'linear','biophys' and 'identity' types of conversion.  The
% linear is a standard convolution.  The biophys is based on Rieke's
% biophysical work.  And identity is a copy operation.
os = osCreate('linear');
% os = osCreate('biophys');

% Compute the photocurrent
os = osCompute(os, absorptions);
 
% Plot the photocurrent for a pixel
% Let's JG and BW mess around with various plotting things to check the
% validity.
%
% osPlot(os,'photo current','cone position',[r,c])
osPlot(os,absorptions);

%% Outer segment: identity for input to RGC
% Input = RGB
osI = osCreate('identity');
osI = osSet(osI, 'rgbData', sceneRGB);
%% Build rgc

% rgc1 = rgcCreate('GLM', scene, sensor, os, 'right', 3.0, 180);

clear params
params.scene = scene; 
params.sensor = absorptions; 
params.outersegment = osI;
params.eyeSide = 'left'; 
params.eyeRadius = 5; 
params.eyeAngle = 90;
rgc1 = rgcCreate('GLM', params);
% rgc1 = rgcCreate('linear', 'sensor', sensor, 'outersegment', os, 'eyeSide','left', 'eyeRadius', 9, 'eyeAngle', 90);

rgc1 = rgcCompute(rgc1, osI);

% rgcPlot(rgc1, 'mosaic');
% rgcPlot(rgc1, 'linearResponse');
rgcPlot(rgc1, 'spikeResponse');