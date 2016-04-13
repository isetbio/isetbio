function iStim = ieStimulusGratingSubunit(varargin)
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
addParameter(p,'barWidth',       3,     @isnumeric);
addParameter(p,'meanLuminance',  200,   @isnumeric);
addParameter(p,'nSteps',         40,    @isnumeric);
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,    @isnumeric);  
addParameter(p,'expTime',        0.01, @isnumeric);
addParameter(p,'timeInterval',   0.01, @isnumeric);

p.parse(varargin{:});
params = p.Results;
fov = params.fov;

wFlag = ieSessionGet('wait bar');

%% Compute a Gabor patch scene as a placeholder for the bar image

% Create display
display = displayCreate('CRT-Sony-HorwitzLab');

% Set up scene, oi and sensor
scene = sceneCreate();
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');
absorptions = sensorCreate('human');
absorptions = sensorSetSizeToFOV(absorptions, fov, scene, oi);

absorptions = sensorSet(absorptions, 'exp time', params.expTime); 
absorptions = sensorSet(absorptions, 'time interval', params.timeInterval); 

%% Compute a dynamic set of cone absorptions for moving bar

fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
% The number of steps must be smaller than the width of the scene
nSteps = params.nSteps;

for t = 1 : nSteps
    if wFlag, waitbar(t/nSteps,wbar); end
        
    sceneSize = sceneGet(scene, 'size');
    barMovie = ones([sceneSize(1)+2*params.barWidth,sceneSize(2), 3])*0.25;  % Gray background
    
    nStripes = floor((sceneSize(1)+params.barWidth)/params.barWidth);
%     randWalk = (round(.25*params.barWidth*randn(1)));
    if t > 1
        % randWalk(t) = sum(randWalk(1:t-1)) + round(.25*params.barWidth);
        randWalk(t) = sum(randWalk(1:t-1)) + (round(.25*params.barWidth*randn(1)));
    else
        randWalk(t) = round(.25*params.barWidth);
    end
    
%     if randWalk(t) > params.barWidth
%         randWalk(t) = params.barWidth;
%     elseif randWalk(t) < -params.barWidth + 1
%         randWalk(t) = -params.barWidth+1;
%     end
    randWalk(t)=t;
    for stripeInd = 1:2:nStripes+2
        % barMovie(randWalk(t)+(stripeInd-1)*params.barWidth+params.barWidth:randWalk(t)+(stripeInd)*params.barWidth+params.barWidth-1,:,:) = 1;          % White bar
        barMovie((stripeInd-1)*params.barWidth+params.barWidth:(stripeInd)*params.barWidth+params.barWidth-1,:,:) = .75;          % White bar
    end
    barMovie = circshift(barMovie,randWalk(t));
    barMovieResize = barMovie(params.barWidth+1:params.barWidth+sceneSize(1),1:sceneSize(2),:);
    
    % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(barMovieResize, 'rgb', params.meanLuminance, display);

    scene = sceneSet(scene, 'h fov', fov);
    if t ==1
        sceneRGB = zeros([sceneGet(scene, 'size'), nSteps, 3]);
    end
    
    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions
    absorptions = sensorCompute(absorptions, oi);

    if t == 1
        volts = zeros([sensorGet(absorptions, 'size') params.nSteps]);
    end
    
    volts(:,:,t) = sensorGet(absorptions, 'volts');
    
    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
absorptions = sensorSet(absorptions, 'volts', volts);
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
