%% t_osIntroduction
%
% In which our heroes lay out the basic architecture of the os class and
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
params.nSteps = 30;     % Number of stimulus frames
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
    % sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
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

%% Show the movie of isomerizations
% 
% % Can we make that movie when we color the cones by type
% 
% vcNewGraphWin;axis image; colormap(gray)
% for ii=1:params.nSteps
%     imagesc(isomerizations(:,:,ii)); pause(.2); 
% end
% 
% % Time series at a point
% vcNewGraphWin; plot(squeeze(isomerizations(1,1,:)))

%% Movie of the cone absorptions over cone mosaic
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

%% Outer segment calculation for the linear model

% The outer segment converts cone absorptions into cone photocurrent.
% There are 'linear','biophys' and 'identity' types of conversion.  The
% linear is a standard convolution.  The biophys is based on Rieke's
% biophysical work.  And identity is a copy operation.
os = osCreate('linear');
 
% Compute the photocurrent
os = osCompute(os, absorptions);
 
% Plot the photocurrent for a pixel
% Let's JG and BW mess around with various plotting things to check the
% validity.
%
% osPlot(os,'photo current','cone position',[r,c])
osPlot(os,absorptions);

%% Examine isomerizations with movies and plots
% Show an image of the cone mosaic
cone_mosaic = sensorGet(absorptions,'cone type');

figure; 
imagesc(5-cone_mosaic);
colormap(jet(4));
title('cone mosaic','fontsize',14)

% Show isomerizations movie
% Note black pixels which are S cones that do not absorb many photons
% in comparison to L and M cones
isomerizations = sensorGet(absorptions,'photons');
figure;
for frame1 = 1:30
    
    imagesc(isomerizations(:,:,frame1));
    colormap gray; 
    drawnow;
end
title('t\_osIntroduction isomerizations');

% Plot isomerizations over time
% Find coordinates of L, M and S cones, get photon signals.
figure; cind = 'rgb';
for cone_type = 2:4
    [cone_locations_x cone_locations_y] = find(cone_mosaic==cone_type);
    
    for cone_number = 1:10%length(cone_locations_x)
        % figure;
        hold on;
        plot(2*rand+squeeze(isomerizations(cone_locations_x(cone_number), cone_locations_y(cone_number),:)),cind(cone_type-1));
    end
end
title('L, M and S cone photons');
 xlabel('Frame number'); ylabel('Photons');
%% Examine linear response with movies and plots
% Show photocurrents movie
coneCurrentSignal = osGet(os,'coneCurrentSignal');
figure;
% cone_type = 4;
% clear cone_locations_x cone_locations_y
% [cone_locations_x cone_locations_y] = find(cone_mosaic==cone_type);
for frame1 = 1:30
   
    imagesc(coneCurrentSignal(:,:,frame1));
%     coneFrame = zeros(size(coneCurrentSignal,1),size(coneCurrentSignal,2));
%     coneFrame(cone_locations_x, cone_locations_y,frame1) = coneCurrentSignal(cone_locations_x,cone_locations_y,frame1);
%     imagesc(coneFrame(:,:,frame1));
    colormap gray; 
%     caxis([-70 -40]);
    drawnow;
end
title('t\_osIntroduction linear photocurrent');
% close;


% Plot photocurrents over time
% Find coordinates of L, M and S cones, get current signals.
figure; cind = 'rgb';
for cone_type = 2:4
[cone_locations_x cone_locations_y] = find(cone_mosaic==cone_type);

for cone_number = 1:10%length(cone_locations_x)
    % figure;
    hold on;
    plot(2*rand+squeeze(coneCurrentSignal(cone_locations_x(cone_number), cone_locations_y(cone_number),:)),cind(cone_type-1));
end
end
xlabel('Frame number'); ylabel('Current (pA)');
title('L, M and S cone currents');

%% Rieke biophysics case

% Create the outer segment structure
osB = osCreate('biophys');
 
% Compute the photocurrent from the absorptions
osB = osCompute(osB, absorptions);
 
% Plot the photocurrent for a pixel
% Let's JG and BW mess around with various plotting things to check the
% validity.

osPlot(osB,absorptions)

%% Examine biophysical response with a movie

coneCurrentSignalB = osGet(osB,'coneCurrentSignal');
figure;
cone_type = 2;
for frame1 = 1:30
   
    imagesc(coneCurrentSignalB(:,:,frame1));
    colormap gray;
%     caxis([-70 0]);
    drawnow;
end
% close;
title('t\_osIntroduction biophys photocurrent');