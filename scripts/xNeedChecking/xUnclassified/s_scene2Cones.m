% s_rgcScene2ConesEyeMovements
%
% The calculation from a scene to the quantum absorptions in a cone array,
% allowing for small, fixational eye movements
%
% See also: s_rgcScene2Cones (without eye movments) and s_rgcCones2RGC for
% the next steps.
%
% (c) Stanford VISTA Team

%%
s_initISET;


%% User-adjustable parameters
% The script is getting long. Let's list all adjudtable parameters at the
% top. The rest of the script will do the computations and display the
% resutls.

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
% We assume the eye makes small movements around fixation.


% --- Cone Adapation ---
% Now that we simulate eye movmements, perhaps we shoud do adaptation by
% individual cone. If there there are no eye movements, then we would like
% image to fade (cone image becomes uniform, other than noise)
typeAdapt = 'rieke';

%% SCENE
% Create scene from image
scene = sceneFromFile(im, 'rgb', meanL, dCal);
scene = sceneSet(scene,'fov',fov);      % adjust fov if specified
scene = sceneSet(scene,'distance',vd);  % Two meters
if ~isempty(illuminant)
    scene = sceneAdjustIlluminant(scene, illuminant);
end

vcAddAndSelectObject(scene); sceneWindow;

%% Human Optics
% Create human optics
oi = oiCreate('wvf human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;


%% Human Sensor
% initialize cone parameters
coneP = coneCreate;
coneP.rgbDensities = coneDensity;
coneP.coneAperture = coneAperture;

% Create the sensor
sensor = sensorCreate('human', coneP);
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

% TODO: Set exposure time and eye movement

% Set voltage swing
sensor = sensorSet(sensor,'pixel voltage swing', voltageSwing);

%% Cone Absorptions
% The coneAbsorptions routine handles eye movements also.  If eyemovement
% parameters were set, the returned sensor will be a row x col x nFrames
% matrix. Each frame corresponds to the absorptions in one movement
% position.
sensor = coneAbsorptions(sensor, oi);
volts  = sensorGet(sensor,'volts');

%%  Adaptation
% The cone output voltages normally operate in a specific range specified
% by the voltage swing of the receptor.  Suppose we set the voltage range
% for a cone to 0-1V.  Then as the scene illumination level swings across
% many orders of magnitude, we still need the cones to be operating with
% some contrast inside of their 0-1V range.  Otherwise the RGC won't see
% the signals.
%
% In general, we need a model of the cone voltage output that is better
% than the linear model followed by a gain control (adpatation) that we use
% here.  It could be a log-model, for example.
%
% See coneAdapt for details
[~, adaptedCur] = coneAdapt(sensor, typeAdapt);

%% End here

% Debugging below

%%  Examine the data using s_rgcCones2RGC

% We should create cone visualization routines that don't require going to
% the rgc structure.
% Make a coneVisualize() routine and separate stuff out?

% creating the RGC parameter object
rgcP = rgcParameters;

% The absorptions structure is also the RGC parameter data
% rgcP.set('cone voltages',cones.data);
% rgcP.set('sensor',cones.sensor);
% rgcP.set('oi',cones.oi);
rgcP.set('cone voltages',adaptedCur); % this should be cone current
rgcP.set('sensor',sensor);
rgcP.set('oi',oi);


% What is the default?  RGC spacing the same as cone spacing (I think).
rgcP.addLayer('on parasol', 20);  

rgcComputeSpikes(rgcP);

%% Visualize 

% NEED TO USE implay rather than mplay.  Re-write required.

rgcVisualize('Cone Voltages',rgcP); 
rgcVisualize('Linear Timeseries',rgcP); 

%% End