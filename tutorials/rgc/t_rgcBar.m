%% t_rgcBar
% 
% A tutorial for the isetbio RGC object. A moving bar is created 
% in isetbio as a scene, and the sensor, outer segment and RGC objects are
% computed in response to the scene. 
% 
% (JRG) modified 10/2015

%% Init
% ieInit;

%% Compute a Gabor patch scene

% Set up Gabor stimulus using sceneCreate('harmonic',params)
fov = 1.2;

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

params.expTime = 0.001;
params.timeInterval = 0.001;
params.nSteps = 200;     % Number of stimulus frames


%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');

sensor = sensorCreate('human');
% 
% % paramsStim.fov = 8;
% paramsStim.radius = 3.5;
% paramsStim.theta = 330;
% paramsStim.side = 'left';
% 
% if paramsStim.radius == 0
%     
%     sensor = sensorCreate('human');
% 
% else
%     
%     coneP = coneCreate; % The cone properties properties
%     retinalRadiusDegrees = paramsStim.radius;
%     retinalPolarDegrees  = paramsStim.theta;
%     whichEye             = paramsStim.side;
%     retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; 
%     sensor = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
% end

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
    barMovie(:,1+floor(t/1):floor(t/1)+params.barwidth,:) = 0.5 + 0.499*ones(params.row,params.barwidth,3);
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

%% Outer segment calculation - linear

% Input = RGB
% osL = osCreate('linear');
osL = osCreate('biophys');

% Set cone spacing so that 
arrayWidth = sensorGet(sensor,'width');
osL = osSet(osL, 'patch size', arrayWidth);

% The size of the whole mosaic
% I think there is a function that gets this value, maybe in oiGet.
%  umPerDeg = 300;
%  coneSpacing = scene.wAngular*umPerDeg;

tSampling = sensorGet(sensor,'time interval','sec');
osL = osSet(osL, 'time step', tSampling);

clear paramsOSL
% Set circular convolution, only steady state
% paramsOSL.convolutionType = 1; 
paramsOSL.ecc = 20; % mm
paramsOSL.singleType = 1;

sensorVolts = sensorGet(sensor,'volts');
paramsOSL.bgVolts = 1*mean(sensorVolts(:));

% osL = osSet(osL, 'noiseflag',1);

% os = osSet(os, 'rgbData', sceneRGB);
osL = osCompute(osL, sensor, paramsOSL);

% Plot the photocurrent for a pixel
osPlot(osL,sensor);

%% Find bipolar responses
clear bp 

bp = bipolar(osL,'filterType',1);
% bp = bipolar(os,cellType,2);

bp = bipolarCompute(bp, osL);

% bipolarPlot(bp,'response');
%% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaSU = irCreate(bp, params);

% innerRetinaSU.mosaicCreate('model','glm','type','off midget');
innerRetinaSU.mosaicCreate('model','glm','type','on midget');
% innerRetinaSU.mosaicCreate('model','glm','type','off parasol');
% innerRetinaSU.mosaicCreate('model','glm','type','on midget');
%% Compute RGC response

innerRetinaSU = irCompute(innerRetinaSU, bp);

%%
% irPlot(innerRetinaSU, 'mosaic');
% irPlot(innerRetinaSU, 'linear');
% irPlot(innerRetinaSU, 'raster');
irPlot(innerRetinaSU, 'psth','type',1);
%%
for i1 = 1:10; 
    for j1 = 1:10; 
        hold on;
        subplot(10,10,(j1-1)*10+i1); 
        hold on;
        axis([0 1.1 0 50]); 
    end; 
end
%% Outer segment calculation
% 
% Input = RGB
os = osCreate('displayrgb');

coneSpacing = sensorGet(sensor,'width','um');
coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', coneSpacing);

coneSampling = sensorGet(sensor,'time interval','sec');
os = osSet(os, 'timeStep', .008);

os = osSet(os, 'rgbData', sceneRGB);
% os = osCompute(sensor);

% % Plot the photocurrent for a pixel
% osPlot(os,sensor);

%% Outer segment calculation - biophysical

% os = osCreate('biophys');
% 
% % Set cone spacing so that 
% arrayWidth = sensorGet(sensor,'width');
% os = osSet(os, 'patch size', arrayWidth);
% 
% % The size of the whole mosaic
% % I think there is a function that gets this value, maybe in oiGet.
% %  umPerDeg = 300;
% %  coneSpacing = scene.wAngular*umPerDeg;
% 
% % Change this to .001 for stability
% tSampling = sensorGet(sensor,'time interval','sec');
% os = osSet(os, 'time step', tSampling);
% 
% % os = osSet(os, 'rgbData', sceneRGB);
% os = osCompute(os, sensor);
% 
% % Plot the photocurrent for a pixel
% osPlot(os,sensor);
%% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina = irCreate(os, params);

innerRetina.mosaicCreate('model','glm','type','on midget');
%% Compute RGC response

innerRetina = irCompute(innerRetina, os);

%%
% irPlot(innerRetina, 'mosaic');
% irPlot(innerRetina, 'linear');
% irPlot(innerRetina, 'raster');
irPlot(innerRetina, 'psth');
