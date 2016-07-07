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
%   clear params; params.barWidth = 10; params.fov=0.6;
%   iStim = ieStimulusBar(params);
%   coneImageActivity(iStim.absorptions,'dFlag',true);
%   vcAddObject(iStim.absorptions); sensorWindow;
% 
% 3/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser;
addParameter(p,'barWidth',       5,     @isnumeric);
addParameter(p,'meanLuminance',  200,   @isnumeric);
addParameter(p,'nSteps',         [],    @isnumeric);
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,    @isnumeric);  
addParameter(p,'expTime',        0.005, @isnumeric);
addParameter(p,'timeInterval',   0.005, @isnumeric);
addParameter(p,'display',   'LCD-Apple',@ischar);

% Retinal patch parameters
addParameter(p,'radius',            0,  @isnumeric);
addParameter(p,'theta',            0,  @isnumeric);
addParameter(p,'side',            'left',  @ischar);

p.parse(varargin{:});
params = p.Results;
fov = params.fov;

wFlag = ieSessionGet('wait bar');

%% Compute a Gabor patch scene as a placeholder for the bar image

% Create display
display = displayCreate(params.display);

% Set up scene, oi and sensor
scene = sceneCreate();
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');
% absorptions = sensorCreate('human');
% absorptions = sensorSetSizeToFOV(absorptions, fov, scene, oi);
% 
% absorptions = sensorSet(absorptions, 'exp time', params.expTime); 
% absorptions = sensorSet(absorptions, 'time interval', params.timeInterval); 


if params.radius == 0
    
    absorptions = sensorCreate('human');

else
    
    coneP = coneCreate; % The cone properties properties
    retinalRadiusDegrees = params.radius;
    retinalPolarDegrees  = params.theta;
    whichEye             = params.side;
    retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; 
    absorptions = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
end

absorptions = sensorSetSizeToFOV(absorptions, params.fov, scene, oi);

% Set aspect ratio
% At this point, we have the right cone density and the right number in
% cols, now we just need to set the rows to have the same aspect ratio as
% the input movie.
% sensorSize = sensorGet(absorptions,'size');
% aspectRatioMovie = size(movieInput,1)/size(movieInput,2);
% absorptions = sensorSet(absorptions,'size',[aspectRatioMovie*sensorSize(2) sensorSize(2)]);


%% Compute a dynamic set of cone absorptions for moving bar

fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
% The number of steps must be smaller than the width of the scene
nSteps = params.nSteps;
if isempty(nSteps), 
    nSteps = sceneGet(scene,'cols') - params.barWidth; 
end
nSteps = min(sceneGet(scene,'cols') - params.barWidth, nSteps);

cMosaic = coneMosaic;

for t = 1 : nSteps
    if wFlag, waitbar(t/nSteps,wbar); end
        
    barMovie = ones([sceneGet(scene, 'size'), 3])*0.001;  % Gray background
    barMovie(:,t:(t+params.barWidth-1),:) = 1;          % White bar

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
    
    % Compute absorptions
    % absorptions = sensorCompute(absorptions, oi);
    cMosaic.compute(oi);
    
    if t == 1
        absorptions = zeros([cMosaic.mosaicSize, params.nSteps]);
    end
    
    absorptions(:,:,t) = cMosaic.absorptions;
    
    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
% absorptions = sensorSet(absorptions, 'volts', absorptions);
% vcAddObject(sensor); sensorWindow;

% These are both the results and the objects needed to recreate this
% script. So calling isomerizationBar(iStim) should produce the same
% results.

iStim.params  = params;
iStim.display = display;
iStim.scene   = scene;
iStim.sceneRGB = sceneRGB;
iStim.oi      = oi;
iStim.absorptions  = absorptions;
end
