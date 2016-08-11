function iStim = ieStimulusBar(varargin)
% Creates a movie/dynamic scene stimulus in isetbio of a white bar on a
% black background that sweeps from left to right.
% 
% Inputs: a structure that defines the parameters of the bar stimulus.
% 
% Outputs: iStim is a structure that contains the display, the scene, the
%   optical image and the sensor.
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

p = inputParser;
addParameter(p,'barWidth',       5,     @isnumeric);
addParameter(p,'meanLuminance',  200,   @isnumeric);
addParameter(p,'nSteps',         inf,   @isnumeric); % determined by cols
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,   @isnumeric);  
addParameter(p,'expTime',        0.005, @isnumeric);
addParameter(p,'timeInterval',   0.005, @isnumeric);
addParameter(p,'display',   'LCD-Apple',@ischar);

% Retinal patch parameters
addParameter(p,'radius',         0,  @isnumeric);
addParameter(p,'theta',          0,  @isnumeric);
addParameter(p,'side',           'left',  @ischar);

p.parse(varargin{:});
params = p.Results;
fov = params.fov;

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

cm = coneMosaic;
cm.pigment.width = coneSz(1); cm.pigment.height = coneSz(2);

% set size to field of view
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cm.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

%% Compute a dynamic set of cone absorptions for moving bar
fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
wbar = waitbar(0,'Stimulus movie');
% Loop through frames to build movie
% The number of steps must be smaller than the width of the scene
nSteps = min(sceneGet(scene,'cols') - params.barWidth, params.nSteps);
for t = 1 : nSteps
    waitbar(t/nSteps,wbar);
        
    barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;  % Gray background
    barMovie(:,t:(t+params.barWidth-1),:) = 1;            % White bar

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
