function iStim = ieStimulusGaborColor(varargin)
% Creates a movie/dynamic of the cone absorptions of a drifting Gabor patch 
% 
%    ieStimulusGabor(varargin)
%
% The Gabor image is created with an equal photon spectral power
% distribution.  The parameters of the Gabor are set by the input structure
% (or key/value parameters) below.
%
% Inputs: 
%   pGabor: parameter structure that defines the Gabor stimulus parameters
%           The parameters and their defaults are 
%
% 'row',              4 x Nyquist   image size
% 'col',              4 x Nyquist   We multiply 4 times fov times freq
% 'meanLuminance',  200   (cd/m2)
% 'freq',            6    Spatial frequency c/deg
% 'contrast'         1    Harmonic contrast
% 'ph',              0    Phase of the harmonic
% 'ang',             0    Angle (0 = horizontal) of the grating variation
% 'GaborFlag',       1    Std of Gaussian where 1 means min(row,col) 
% 'fov',            0.6   Field of view
% 'expTime',        0.005 Exposure time of the sensor
% 'nCycles',        4     Number of cycles through the harmonic
% 'nSteps',         15*4    Total number of steps for all the harmonic cycles
%                         There are nSteps/nCycles per each cycle.
% 'distance'        0.3   Meters from the viewing screen
%
% Outputs: 
%   iStim - a structure that contains
%     params used to create this
%     scene (first frame)
%     sceneRGB
%     optical image (first frame)
%     human cone absorptions (dynamic)
% 
% Examples:
% Coarsely stepped, pretty tight Gabor window
%   nSteps = 20; GaborFlag = 0.2; fov = .5;
%   iStim = ieStimulusGabor('nSteps',nSteps,'GaborFlag',GaborFlag,'fov',fov);
%   m = coneImageActivity(iStim.absorptions,'dFlag',true);
%   sceneShowImage(iStim.scene);
%
%  Higher spatial frequency, more steps
%   params.freq = 6; params.nSteps = 50; params.GaborFlag = 0.2;
%   iStim = ieStimulusGabor(params);
%   dFlag.vname = 'Gabor_6f';
%   dFlag.FrameRate = 30;
%   coneImageActivity(iStim.absorptions,'dFlag',dFlag);
%   sceneShowImage(iStim.scene);
%
% 3/2016 JRG (c) isetbio team

%% Parse inputs
p = inputParser;

addParameter(p,'meanLuminance',  100,   @isnumeric);
addParameter(p,'nSteps',         60,    @isnumeric);
addParameter(p,'row',             0,    @isnumeric);  
addParameter(p,'col',             0,    @isnumeric);  
addParameter(p,'freq',            6,    @isnumeric);  
addParameter(p,'contrast',        1,    @isnumeric);  
addParameter(p,'ph',              0,    @isnumeric);  
addParameter(p,'ang',             0,    @isnumeric);  
addParameter(p,'GaborFlag',       1,    @isnumeric);  
addParameter(p,'expTime',        0.005, @isnumeric);
addParameter(p,'nCycles',        4, @isnumeric);
addParameter(p, 'colorInd',      2, @isnumeric);

% Viewing parameters
addParameter(p,'fov',            0.6, @isnumeric);
addParameter(p,'distance',       0.3, @isnumeric);   % Distance to screen

% Field of view
p.parse(varargin{:});
paramsIn = p.Results;

fov = paramsIn.fov;
if paramsIn.row == 0
    % Make sure we have enough row and column samples to avoid aliasing the
    % frequency.  Four is arbitrary, but twice Nyquist.
   paramsIn.row = 4*paramsIn.freq*paramsIn.fov;
   paramsIn.col = 4*paramsIn.freq*paramsIn.fov;
end

%%

nContrast   = 4;
maxContrast = linspace(0,1,nContrast);

% contrastArr = [0,1];%
% contrastArr = (0:0.25:1); % must start with 0
% nContrast = length(contrastArr);

noiseIterations = 100;    % more iterations will improve accuracy but take longer!
pooledData      = cell(1,nContrast);
rocArea         = zeros(1,nContrast);

% Load the display from the Horwitz Lab
display = displayCreate('CRT-Sony-HorwitzLab');

% Basic human optics parameters.  Perhaps we need to get the macaque optics
% built here.
oi  = oiCreate('wvf human');


% Set parameters for building the dynamic scene/oi/sensor The
% stimulus consists of a color opponent Gabor patch whose phase
% varies over time. The stimulus parameters specified in Fig. 6 of
% the manuscript.

% parameters found in Fig. 6 caption
params = paramsGaborColorOpponent();
params.contrast = 0;                  % set max contrast of gabor
params.color = paramsIn.colorInd;    % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
params.image_size = 64;               % scene is (image_size X image_size) pixels
params.fov        = 0.6;
params.nSteps     = 60; % 666; % 160+346+160
params.contrast = paramsIn.contrast;%maxContrast(contrastInd);

% This is the contrast at each step in the simulation.  It ramps up
% to the max (1 in this case), stays there, and then ramps back
% down.
%         timeContrast = ...
%             [linspace(0,1,(1/4)*params.nSteps), ...
%              ones(1,params.nSteps/2),...
%              linspace(1,0,(1/4)*params.nSteps)];
timeContrast = ones(params.nSteps,1);
% plot(timeContrast)

% We build a dummy scene here just so we can subsequently calculate
% the sensor size.  But this scene itself is not used.  Rather we
% build a series of scenes below.
stimulusRGBdata = imageHarmonicColor(params); % sceneCreateGabor(params);
scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor

coneP = coneCreate; % The cone properties properties

% see caption for Fig. 4 of Horwitz, Hass, Rieke, 2015, J. Neuro.
%         retinalPosDegAz = 5; retinalPosDegEl = -3.5;
%         retinalRadiusDegrees = sqrt(retinalPosDegAz^2+retinalPosDegEl^2);
%         retinalPolarDegrees  = abs(atand(retinalPosDegEl/retinalPosDegAz));
%         retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; whichEye = 'right';
%         sensor = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, params.fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', params.expTime); % 1 ms
sensor = sensorSet(sensor, 'time interval', params.timeInterval); % 1 ms

% This computes with no sensor or photon noise, just the mean
sensor = sensorSet(sensor,'noise flag',0);

% According to the paper, cone collecting area is 0.6 um^2
% sensor = sensorSet(sensor 'pixel pd width', 0.774e-6); % photo-detector width
% sensor = sensorSet(sensor, 'pixel pd height', 0.774e-6);

% Macular pigment is attached to the sensor
% macular pigment absorbance was scaled to 0.35 at 460 nm
% macular = sensorGet(sensor, 'human macular');
% macular = macularSet(macular, 'density', 0.35);
% sensor = sensorSet(sensor, 'human macular', macular);

% Compute a dynamic set of cone absorptions
% ieSessionSet('wait bar',true);
fprintf('Computing dynamic scene/oi/sensor data    ');

% if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
for t = 1 : params.nSteps
    
    fprintf('\b\b\b%02d%%', round(100*t/params.nSteps));
    
%     if wFlag, waitbar(t/params.nSteps,wbar); end
    
    % Update the phase of the Gabor
    params.ph = 2*pi*(t-1)/params.nSteps; % one period over nSteps
    params.contrast = params.contrast*timeContrast(t);
    
    % scene = sceneCreate('harmonic', params);
    % scene = sceneSet(scene, 'h fov', fov);
    
    % The imageHarmonicColor should be moved to ISETBIO.
    stimulusRGBdata = imageHarmonicColor(params); % sceneCreateGabor(params);
    
    scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
    scene = sceneSet(scene, 'h fov', params.fov);
    
    % The mean scene luminance is set to 200 cd/m2, based on the
    % calibration of the monitor
    scene = sceneAdjustLuminance(scene, 200);
    
    % oi  = oiCreate('wvf human');
    % Compute optical image
    oi = oiCompute(oi, scene);
    
    % Compute absorptions
    sensor = sensorCompute(sensor, oi);
    
    if t == 1
        volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
        
    end
    
    volts(:,:,t) = sensorGet(sensor, 'volts');
    
    % In case we need to visualize, you can do this on each
    % iteration
    % vcAddObject(scene); sceneWindow
    % pause(.1);
end
fprintf('\n');
toc
% if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
sensor = sensorSet(sensor, 'volts', volts);


% Save all the inputs to rerun 
iStim.params   = params;     % Parameters to rerun this function
iStim.scene    = scene;      % Base scene
% iStim.sceneRGB = sceneRGB;   % Used for identity case.
iStim.oi       = oi;         % 
iStim.sensor   = sensor;
end