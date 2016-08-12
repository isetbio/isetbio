function iStim = ieStimulusBar(varargin)
% Creates a dynamic cone mosaic response to a moving bar stimulus
% 
% Inputs: a parameter structure that defines
%   * the bar stimulus properties
%   * cone mosaic properties
%   
%  
%  Stimulus parameters:
%    display, barWidth, meanLuminance, nSteps,row, col, fov, grayStart,
%    grayEnd
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
%   iStim.cMosaic.window;
%  Returns the same 
%   foo = ieStimulusBar(iStim.params);
%   foo.cMosaic.window;
%
% 3/2016 JRG (c) isetbio team

%% Parse inputs

% Loop through frames to build movie
% The number of steps must be smaller than the width of the scene
grayStart = 50; grayEnd = 20;


p = inputParser;
% Stimulus parameters
addParameter(p,'display',   'LCD-Apple',@ischar);
addParameter(p,'barWidth',       5,     @isnumeric); % Pixels
addParameter(p,'meanLuminance',  200,   @isnumeric); % Cd/m2
addParameter(p,'nSteps',         inf,   @isnumeric); % determined by cols
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,   @isnumeric); % Deg 
addParameter(p,'grayStart',      75,    @isnumeric); % ms 
addParameter(p,'grayEnd',        50,    @isnumeric); % ms 

% OS and mosaic parameters
addParameter(p,'os',            'linear',@ischar);
addParameter(p,'expTime',        0.005, @isnumeric); % Sec
addParameter(p,'radius',         0,  @isnumeric);  % Degrees?
addParameter(p,'theta',          0,  @isnumeric);  % Degrees?
addParameter(p,'side',           'left',  @ischar);% Left/right

p.parse(varargin{:});
params = p.Results;
fov = params.fov;
osType = p.Results.os;

% Turn off so we supercede with this wait bar
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
coneD = coneDensity(eccMM, [params.radius params.theta], params.side);
coneSz(1) = sqrt(1./coneD) * 1e-3;  % avg cone size with gap in meters
coneSz(2) = coneSz(1);

if strcmpi(osType, 'biophys');
    osCM = osBioPhys();
    cm = coneMosaic('os',osCM);
else
    cm = coneMosaic;
end

cm.pigment.width = coneSz(1); cm.pigment.height = coneSz(2);

% set size to field of view
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cm.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

%% Compute a dynamic set of cone absorptions for moving bar
fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
wbar = waitbar(0,'Stimulus movie');

% nSteps = min(sceneGet(scene,'cols')+grayStart+grayEnd, params.nSteps);
nStepsStim = (sceneGet(scene,'cols')+grayStart-params.barWidth);
nSteps = nStepsStim + grayEnd;
for t = 1 : nSteps
    waitbar(t/nSteps,wbar);
        
    barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;  % Gray background
    
    if t > grayStart && t < nStepsStim
    barMovie(:,t-grayStart+1:(t-grayStart+1+params.barWidth-1),:) = 1;            % White bar
    end
    
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
    
    % Compute absorptions and photocurrent
    cm.compute(oi, 'append', true, 'emPath', [0 0]);
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
iStim.cMosaic  = cm;
end
