function iStim = ieStimulusBar(varargin)
% Creates a dynamic cone mosaic response to a moving bar stimulus
% 
% Inputs: a parameter structure that defines
%   * the bar stimulus properties
%   * cone mosaic properties
%   
%  Stimulus parameters:
%    display, barWidth, meanLuminance, nSteps,row, col, fov, 
%    startFrames   - Frames before stimulus (uniform, mean luminance)
%    stimFrames    - Number of stimulus frames
%    endFrames     - Frames after stimulus end (uniform, mean luminance)
%
%  Cone mosaic parameters:
%    os, expTime, eccentricity, angle, side
%
% Outputs: iStim is a structure that contains 
%   * display model
%   * original static bar scene
%   * optical image of the scene
%   * cone mosaic responses to the scene as the bar translates 
%      (no eye movements)
% 
% Example:
%   clear params; params.barWidth = 10; params.fov=1;
%   iStim = ieStimulusBar(params);
%   iStim.cm.window;
%
%  Returns the same 
%   iStim = ieStimulusBar(iStim.params);
%   iStim.cm.window;

% 3/2016 JRG (c) isetbio team
%
% 08/16/17  dhb  Fix call to coneDensity -> coneDensityReadData.

%% Parse inputs
p = inputParser;

% Stimulus parameters
addParameter(p,'display',   'LCD-Apple',@ischar);
addParameter(p,'barWidth',       5,     @isnumeric); % Pixels
addParameter(p,'meanLuminance',  200,   @isnumeric); % Cd/m2
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,   @isnumeric); % Deg 
addParameter(p,'startFrames',    120,    @isnumeric); % ms 
addParameter(p,'stimFrames',     inf,   @isnumeric); % determined by cols
addParameter(p,'endFrames',      60,    @isnumeric); % ms 
addParameter(p,'cmNoiseFlag',      'none',    @ischar); % 'none','random','frozen' 
addParameter(p,'osNoiseFlag',      'none',    @ischar); % 'none','random','frozen'  
addParameter(p,'integrationTime',    .001,    @isnumeric); % ms 

% OS and mosaic parameters
addParameter(p,'os',            'linear',@ischar);
addParameter(p,'radius',         0,  @isnumeric);  % Degrees?
addParameter(p,'theta',          0,  @isnumeric);  % Degrees?
addParameter(p,'side',           'left',  @ischar);% Left/right

p.parse(varargin{:});
params = p.Results;
fov = params.fov;
osType = p.Results.os;
startFrames = params.startFrames;
endFrames   = params.endFrames;
integrationTime = p.Results.integrationTime;
cmNoiseFlag = p.Results.cmNoiseFlag;
osNoiseFlag = p.Results.osNoiseFlag;

conerows = p.Results.row;
conecols = p.Results.col;

% We insist on turning off the wait bar
wFlag = ieSessionGet('wait bar');
ieSessionSet('wait bar',false);
%% Compute a Gabor patch scene as a placeholder for the bar image

% Create display
display = displayCreate(params.display);

% Set up scene, oi and sensor
scene = sceneCreate();
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

% Set the cone aperture size
cm.pigment.width  = coneSz(1); 
cm.pigment.height = coneSz(2);

% Set cone mosaic field of view to match the scene
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cm.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

% Set the noise flag for the absorptions
cm.noiseFlag = cmNoiseFlag;

%% Compute a dynamic set of cone absorptions for moving bar
fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
wbar = waitbar(0,'Stimulus movie');

% This is the mean oi that we use at the start and end
barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;  % Gray background
scene    = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);
oiMean   = oiCompute(oi,scene);

% nSteps = min(sceneGet(scene,'cols')+grayStart+grayEnd, params.nSteps);
stimFrames = (sceneGet(scene,'cols') - params.barWidth);
nSteps = startFrames + stimFrames + endFrames;
for t = 1 : nSteps
    waitbar(t/nSteps,wbar);

    if ~(t > startFrames && t < (startFrames + stimFrames + 1))
        % Use uniform field oi for time prior to and after stimulus
        oi = oiMean;
    else
        
        % Gray background
        barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;  
        
        % Bar at this time
        colStart = t - startFrames + 1;
        colEnd   = colStart + params.barWidth - 1;
        % barMovie(:,t-startFrames + 1:(t-startFrames+1+params.barWidth-1),:) = 1;
        barMovie(:,colStart:colEnd,:) = 200;

        % Generate scene object from stimulus RGB matrix and display object
        scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);
        
        scene = sceneSet(scene, 'h fov', fov);
        if t ==1
            sceneRGB = zeros([sceneGet(scene, 'size'), nSteps, 3]);
        end
        
        % Get scene RGB data
        sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
        
        % Compute optical image
        oi = oiCompute(oi, scene);
    end   

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

delete(wbar);

% Restore
ieSessionSet('wait bar',wFlag);

% These are both the results and the data needed to run this
% script. So calling isomerizationBar(iStim.params) should produce the same
% results.
iStim.params  = params;
iStim.display = display;
iStim.scene   = scene;
iStim.sceneRGB = sceneRGB;
iStim.oi       = oi;
iStim.cm  = cm;
end
