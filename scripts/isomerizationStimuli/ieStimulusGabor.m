function iStim = ieStimulusGabor(varargin)
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
%   sceneShowImage(iStim.scene);
%   iStim.cMosaic.window;
%
%  Higher spatial frequency, more steps
%   clear params; 
%   params.nCycles = 8;
%   params.freq = 6; params.nSteps = 200; params.GaborFlag = 0.2;
%   iStim = ieStimulusGabor(params);
%   vcAddObject(iStim.oi); oiWindow;
%   iStim.cMosaic.window;

% 3/2016 JRG (c) isetbio team
%
% 08/16/17  dhb  Fix call to coneDensity -> coneDensityReadData.

%% Parse inputs
p = inputParser;

addParameter(p,'display',   'LCD-Apple',@ischar);
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
addParameter(p,'os',            'linear',@ischar);
addParameter(p,'integrationTime',    .001,    @isnumeric); % ms 

addParameter(p,'cmNoiseFlag',      'none',    @ischar); % 'none','random','frozen' 
addParameter(p,'osNoiseFlag',      'none',    @ischar); % 'none','random','frozen' 

% Viewing parameters
addParameter(p,'fov',            0.6, @isnumeric);
addParameter(p,'distance',       0.3, @isnumeric);   % Distance to screen

% Retinal patch parameters
addParameter(p,'radius',         0,  @isnumeric);
addParameter(p,'theta',          0,  @isnumeric);
addParameter(p,'side',           'left',  @ischar);

addParameter(p,'startFrames',    120,    @isnumeric); % ms 
addParameter(p,'stimFrames',     inf,   @isnumeric); % determined by cols
addParameter(p,'endFrames',      60,    @isnumeric); % ms 

% Field of view
p.parse(varargin{:});
params = p.Results;
fov = params.fov;
startFrames = params.startFrames;
endFrames   = params.endFrames;

conerows = p.Results.row;
conecols = p.Results.col;

integrationTime = p.Results.integrationTime;

cmNoiseFlag = p.Results.cmNoiseFlag;
osNoiseFlag = p.Results.osNoiseFlag;

osType = p.Results.os;

if params.row == 0
    % Make sure we have enough row and column samples to avoid aliasing the
    % frequency.  Eight is arbitrary, but four times Nyquist.
   params.row = max(128,8*params.freq*params.fov);
   params.col = max(128,8*params.freq*params.fov);
end

%% Compute a scene

% % Create display
display = displayCreate(params.display);

% Set up scene parameters
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');

% compute cone packing density
fLength = oiGet(oi, 'focal length');
eccMM = 2 * tand(params.radius/2) * fLength * 1e3;
coneD = coneDensityReadData('eccentricity',1e-3*eccMM, 'angle',params.theta, 'whichEye', params.side);
coneSz(1) = sqrt(1./coneD) * 1e-3;  % avg cone size with gap in meters
coneSz(2) = coneSz(1);

if strcmpi(osType, 'biophys');
    osCM = osBioPhys();            % peripheral (fast) cone dynamics
    osCM.set('noise flag','none');
%     osCM = osBioPhys('eccentricity',0);  % foveal (slow) cone dynamics
    cm = coneMosaic('os',osCM);
    
    cm.integrationTime = integrationTime;% cm.os.timeStep;
    cm.os.patchSize = cm.width;
    cm.os.timeStep = cm.integrationTime;
elseif strcmpi(osType,'hex')    
    rng('default'); rng(219347);

    resamplingFactor = 4;
    varyingDensity = false;
    customLambda = 0.6;       % If set to empty, @coneMosaiHex chooses
    cm = coneMosaicHex(resamplingFactor, varyingDensity, customLambda, ...
                             'name', 'the hex mosaic', ...
                             'size', [conerows conecols], ...
                        'noiseFlag', 'none',  ...
                   'spatialDensity', [0 0.62 0.31 0.07] ...
          );
    
    cm.integrationTime = integrationTime;% cm.os.timeStep;
    % Set the mosaic's FOV to a wide aspect ratio
    cm.setSizeToFOVForHexMosaic([0.9 0.6]);
else
    cm = coneMosaic;
    
    cm.integrationTime = integrationTime;% cm.os.timeStep;
end

cm.pigment.width  = coneSz(1);
cm.pigment.height = coneSz(2);

% set size to field of view
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cm.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

% Set the noise flag for the absorptions
cm.noiseFlag = cmNoiseFlag;
%% Compute a dynamic set of cone absorptions
%
% We produce a scene video that translates into an oi video that becomes a
% cone absorption video.
%
% We recreate a set of scenes with different phase positions and produce
% the scenes, ois, and cone absorptions by the loop. The result will be a
% time series of the cone photon absorptions.
%
% We are reluctant to make scene(:,:,:,t) because we are frightened about
% the size.  But it still might be the right thing to do.  So the code here
% is an experiment and we aren't sure how it will go.

wFlag = ieSessionGet('wait bar'); ieSessionSet('wait bar',false);

wbar = waitbar(0,'Stimulus movie');
grayStart = 50; grayEnd = 20;

%%%%%%%%
% This is the mean oi that we use at the start and end
params.ph = 0
scene    = sceneCreate('harmonic', params);
oiMean   = oiCompute(oi,scene);

% nSteps = min(sceneGet(scene,'cols')+grayStart+grayEnd, params.nSteps);
stimFrames = (sceneGet(scene,'cols'));
nSteps = startFrames + stimFrames + endFrames;
%%%%%%%%

% Loop through frames to build movie
for t = 1 : params.nSteps
    waitbar(t/params.nSteps,wbar);
        
    % All we do is update the phase of the Gabor
    params.ph = (2*pi)*params.nCycles*(t-1)/params.nSteps; % one period over nSteps
    
    if t > grayStart && t < params.nSteps-grayEnd
        scene = sceneCreate('harmonic', params);
    else
        contrastOld = params.contrast;
        params.contrast = 0;
        scene = sceneCreate('harmonic', params);
        params.contrast = contrastOld;
    end
    
    scene = sceneAdjustLuminance(scene,params.meanLuminance);
    scene = sceneSet(scene,'distance',params.distance);
    scene = sceneSet(scene, 'h fov', fov);

    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions and photocurrent
    % cm.compute(oi, 'append', true, 'emPath', [0 0]);
    cm.compute(oi, 'emPath', [0 0]);
    if t == 1; absorptionsMat = zeros([size(cm.pattern) nSteps]); end;
    absorptionsMat(:,:,t) = cm.absorptions;
    
end
% Need to compute current after otherwise osAddNoise is wrong

cm.absorptions = absorptionsMat;
cm.emPositions=zeros(nSteps,2);
cm.os.noiseFlag = osNoiseFlag;

if strcmpi(osType, 'biophys');
    
    cm.absorptions = cm.integrationTime*cm.absorptions;
    osBParams.bgR = 10*mean(cm.absorptions(:)./cm.os.timeStep);
    cm.computeCurrent(osBParams);
%     cm.absorptions = cm.integrationTime*cm.absorptions;
%     cm.computeCurrent();
else    
    cm.computeCurrent();
end

delete(wbar); ieSessionSet('wait bar',wFlag);

% Save all the inputs to rerun 
iStim.params   = params;     % Parameters to rerun this function
iStim.scene    = scene;      % Base scene
iStim.oi       = oi;         
iStim.cMosaic  = cm;

end
