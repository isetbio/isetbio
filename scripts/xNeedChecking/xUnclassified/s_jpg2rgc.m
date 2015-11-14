% s_jpg2rgc
%
% The calculation from a scene specified by a jpg to the spiking output of
% RGCs.
% 
%
% See also: t_rgc.m for a basic tutorial on RGC response to a movie
% stimulus.
%
% (c) Stanford VISTA Team

%%
s_initISET;


%% User-adjustable parameters

% ---- Scene ----------
im          = 'eagle.jpg';
fov         =  0.4;             % shrink the fov to speed up calculations
dCal        = 'OLED-Sony';      % display to be used
illuminant  = 'D65.mat';        % illuminant, empty to stick with display
meanL       = 500;              % Make bright for faster run time
vd          = 2;                % Viewing distance in meters

% ---- Sensor ----------
coneAperture   = 2e-6;          % um
coneDensity    = [0 .6 .3 .1];  % proportion of [K,L,M,S]
voltageSwing   = 0.1;           % 100 mV to saturation
readNoise      = 0.0005;        % Noise sd on a read 1 mV

% ---- Eye movements -----
% We assume no eye movements for simplicity at this point.

% --- Cone Adaptation ---

% ---- RGC Response ----

%% SCENE
% Create scene from image
scene = sceneFromFile(im, 'rgb', meanL, dCal);
scene = sceneSet(scene,'fov',fov);      % adjust fov if specified
scene = sceneSet(scene,'distance',vd);  % Two meters
if ~isempty(illuminant)
    scene = sceneAdjustIlluminant(scene, illuminant);
end

% Make sceneRGB a (x,y,t,3) array, where the last dimension is RGB
sceneRGB(:,:,1,:) = sceneGet(scene,'rgb');

vcAddAndSelectObject(scene); sceneWindow;

%% Human Optics
% Create human optics
oi = oiCreate('wvf human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;


%% Human Sensor
% initialize cone parameters
params = coneCreate;
params.rgbDensities = coneDensity;
params.coneAperture = coneAperture;

% Create the sensor
sensor = sensorCreate('human', [], params);
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

% TODO: Set exposure time and eye movement

% Set voltage swing
sensor = sensorSet(sensor,'pixel voltage swing', voltageSwing);

%% Cone Absorptions
% The coneAbsorptions routine handles eye movements also.  If eyemovement
% parameters were set, the returned sensor will be a row x col x nFrames
% matrix. Each frame corresponds to the absorptions in one movement
% position.

% sensor = coneAbsorptions(sensor, oi);
sensor = sensorCompute(sensor, oi);
    
volts  = sensorGet(sensor,'volts');

%% Outer segment - identity

os = osCreate('identity');
os = osSet(os, 'rgbData', sceneRGB);

%% RGC

% Retinal position
eyeSide = 'right';
radiusMicrons = 3;
anglePolarDegrees = 180;

rgc1 = rgcCreate('glm', scene, sensor, os, eyeSide, radiusMicrons, anglePolarDegrees);
rgc1 = rgcSet(rgc1,'numberTrials',10);
rgc1 = rgcCompute(rgc1, os);

rgcPlot(rgc1, 'mosaic');
% rgcPlot(rgc1, 'linearResponse');
% rgcPlot(rgc1, 'spikeResponse');
rgcPlot(rgc1, 'rasterResponse');
% rgcPlot(rgc1, 'psthResponse');


for mosaicInd = 1:5
    
    rasterResponse{mosaicInd} = mosaicGet(rgc1.mosaic{mosaicInd}, 'rasterResponse');
    psthResponse{mosaicInd} = mosaicGet(rgc1.mosaic{mosaicInd}, 'psthResponse');
end

fv = rgcGet(rgc1, 'featureVector');
% The feature vector is currently the average number of spikes per trial
% for each cell. Thus if there are N cells, the feature vector is an
% N-element vector where each entry is the average number of spikes per
% presentation.