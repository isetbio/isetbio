function iStim = ieStimulusGaborColor(varargin)
% Creates a movie/dynamic scene stimulus in isetbio of a white bar on a
% black background that sweeps from left to right.
% 
% Inputs: a structure that defines the parameters of the bar stimulus.
% 
% Outputs: iStim is a structure that contains the display, the scene, the
%   optical image and the sensor.
% 
% Example:
%   clear params; params.barWidth = 10; params.fov=0.6;
%   iStim = ieStimulusBar(params);
%   coneImageActivity(iStim.absorptions,'dFlag',true);
%   vcAddObject(iStim.absorptions); sensorWindow;
% 
% 3/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser;
addParameter(p,'meanLuminance',  100,   @isnumeric);
addParameter(p,'nSteps',         60,   @isnumeric); % determined by cols
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,   @isnumeric);  
addParameter(p,'distance',       0.3,   @isnumeric);   % Distance to screen
addParameter(p,'expTime',        0.005, @isnumeric);
addParameter(p,'timeInterval',   0.005, @isnumeric);
addParameter(p,'display',   'LCD-Apple',@ischar);
addParameter(p,'freq',            6,    @isnumeric);  
addParameter(p,'contrast',        1,    @isnumeric);  
addParameter(p,'ph',              0,    @isnumeric);  
addParameter(p,'ang',             0,    @isnumeric);  
addParameter(p,'GaborFlag',       .15,    @isnumeric);  
addParameter(p,'nCycles',        4,     @isnumeric);
addParameter(p, 'colorInd',      2,     @isnumeric);

% Retinal patch parameters
addParameter(p,'radius',         0,  @isnumeric);
addParameter(p,'theta',          0,  @isnumeric);
addParameter(p,'side',           'left',  @ischar);

p.parse(varargin{:});
params = p.Results;
params.color = params.colorInd;
fov = params.fov;
nSteps = params.nSteps;

if params.row == 0
    % Make sure we have enough row and column samples to avoid aliasing the
    % frequency.  Four is arbitrary, but twice Nyquist.
   params.row = 4*params.freq*params.fov;
   params.col = 4*params.freq*params.fov;
end

%% Generate the display

% Create display
display = displayCreate(params.display);

% Set up scene, oi and sensor
% scene = sceneCreate();
% scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

%% Initialize scene as a Gabor image

% This is the contrast at each step in the simulation.  It ramps up
% to the max (1 in this case), stays there, and then ramps back
% down.
timeContrast = ones(params.nSteps,1);

% We build a dummy scene here just so we can subsequently calculate
% the sensor size.  But this scene itself is not used.  Rather we
% build a series of scenes below.
stimulusRGBdata = imageHarmonicColor(params); % sceneCreateGabor(params);
scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
scene = sceneSet(scene, 'h fov', fov);

% vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');
% absorptions = sensorCreate('human');
% absorptions = sensorSetSizeToFOV(absorptions, fov, scene, oi);
% 
% absorptions = sensorSet(absorptions, 'exp time', params.expTime); 
% absorptions = sensorSet(absorptions, 'time interval', params.timeInterval); 

% compute cone packing density
fLength = oiGet(oi, 'focal length');
eccMM = 2 * tand(params.radius/2) * fLength * 1e3;
coneD = coneDensity(eccMM, [params.radius params.theta], params.side);
coneSz = sqrt(1./coneD) * 1e-3;  % avg cone size with gap in meters

cm = coneMosaic;
cm.pigment.width = coneSz(1); cm.pigment.height = coneSz(2);

% set size to field of view
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cm.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

%% Compute a dynamic set of cone absorptions for moving bar
fprintf('Computing cone isomerizations:    \n');

% Loop through frames to build movie
absorptions = zeros([cm.mosaicSize nSteps]);
for t = 1 : nSteps
        
    % Update the phase of the Gabor
    params.ph = 2*pi*(t-1)/params.nSteps; % one period over nSteps
    params.contrast = params.contrast*timeContrast(t);
        
    stimulusRGBdata = imageHarmonicColor(params); 
    
    % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
    scene = sceneSet(scene, 'h fov', fov);
    
    % The mean scene luminance is set to 200 cd/m2, based on the
    % calibration of the monitor
    scene = sceneAdjustLuminance(scene, 200);
    
    if t ==1
        sceneRGB = zeros([sceneGet(scene, 'size'), nSteps, 3]);
    end
    
    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions
    absorptions(:,:,t) = cm.compute(oi, 'currentFlag', false);
end

% Set the stimuls into the sensor object
cm.emPositions = zeros(nSteps, 2);
cm.absorptions = absorptions;

cm.compute(oi, 'currentFlag', true);

% These are both the results and the objects needed to recreate this
% script. So calling isomerizationBar(iStim) should produce the same
% results.
iStim.params  = params;
iStim.display = display;
iStim.scene   = scene;
iStim.sceneRGB = sceneRGB;
iStim.oi       = oi;
iStim.cMosaic  = cm;
end
