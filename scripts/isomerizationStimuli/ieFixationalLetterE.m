function iStim = ieFixationalLetterE(varargin)
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

%%

% stimSize = 32;
% stimFrames = 40;
% stim = zeros(stimSize,stimSize,stimFrames);
% 
% eVert = 6;
% 
% figure;
% for t = 1:stimFrames
% %     stimFrame = 
%     stim(:,:,t) = zeros(stimSize,stimSize);
%     jitter = round((stimSize/4)*(rand(1,2)-0.5));
%     stim(round(stimSize/2)-eVert+1+jitter(1):round(stimSize/2)+eVert+jitter(1),round(stimSize/2)+jitter(2),t) = ones(size(stimSize-eVert+1:stimSize+eVert));
%     stim(round(stimSize/3)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),t) = ones(size(round(stimSize/2)+1:round(stimSize/2)+eVert));
%     stim(round(stimSize/2)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),t) = ones(size(round(stimSize/2)+1:round(stimSize/2)+eVert));
%     stim(round(stimSize*2/3)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),t) = ones(size(round(stimSize/2)+1:round(stimSize/2)+eVert));
%     imagesc(stim(:,:,t)); drawnow
% end

%% Parse inputs

p = inputParser;
addParameter(p,'eWidth',         2,     @isnumeric);
addParameter(p,'eHeight',        9,     @isnumeric);
addParameter(p,'meanLuminance',  200,   @isnumeric);
addParameter(p,'nSteps',         20,    @isnumeric);
addParameter(p,'row',            64,    @isnumeric);  
addParameter(p,'col',            64,    @isnumeric);  
addParameter(p,'fov',            0.6,    @isnumeric);  
addParameter(p,'expTime',        0.005, @isnumeric);
addParameter(p,'timeInterval',   0.005, @isnumeric);

p.parse(varargin{:});
params = p.Results;
fov = params.fov;
eWidth = params.eWidth;
eHeight = params.eHeight;

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
        
    stimDimensions = sceneGet(scene, 'size'); stimSize = stimDimensions(1);
    % eHeight = 13; eWidth = 3;
    eMovie = zeros(stimSize,stimSize,3);
    jitter = round((stimSize/4)*(rand(1,2)-0.5));
    
    %        eVert = 10;
    %     eMovie = zeros(stimSize,stimSize,3);
    %     jitter = round((stimSize/4)*(rand(1,2)-0.5));
    %     eMovie(round(stimSize/2)-eVert+1+jitter(1):round(stimSize/2)+eVert+jitter(1),round(stimSize/2)+jitter(2),:) = ones([size(stimSize-eVert+1:stimSize+eVert),3]);
    %     eMovie(round(stimSize/3)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),:) = ones([size(round(stimSize/2)+1:round(stimSize/2)+eVert),3]);
    %     eMovie(round(stimSize/2)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),:) = ones([size(round(stimSize/2)+1:round(stimSize/2)+eVert),3]);
    %     eMovie(round(stimSize*2/3)+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eVert+jitter(2),:) = ones([size(round(stimSize/2)+1:round(stimSize/2)+eVert),3]);
    
    eMovie(round(stimSize/2)-eHeight+1+jitter(1):round(stimSize/2)+eHeight+jitter(1),...
        round(stimSize/2)-eWidth+1+jitter(2):round(stimSize/2)+eWidth+jitter(2),:) = ...
        ones([size(stimSize-eHeight+1:stimSize+eHeight),...
        size(round(stimSize/2)-eWidth+1+jitter(2):round(stimSize/2)+eWidth+jitter(2)),3]);
    
    eMovie(round(stimSize/2)-eWidth+1+jitter(1):round(stimSize/2)+eWidth+jitter(1),round(stimSize/2)+1+jitter(2):round(stimSize/2)+eHeight+jitter(2),:) = ones([length(round(stimSize/2)-eWidth+1+jitter(1):round(stimSize/2)+eWidth+jitter(1)),length(round(stimSize/2)+1:round(stimSize/2)+eHeight),3]);
    eMovie(round(stimSize/3)-eWidth+1+jitter(1):round(stimSize/3)+eWidth+jitter(1),round(stimSize/2)-1+jitter(2):round(stimSize/2)+eHeight+jitter(2),:) = ones([length(round(stimSize/3)-eWidth+1+jitter(1):round(stimSize/3)+eWidth+jitter(1)),length(round(stimSize/2)-1:round(stimSize/2)+eHeight),3]);
    eMovie(round(2*stimSize/3)-eWidth+1+jitter(1):round(2*stimSize/3)+eWidth+jitter(1),round(stimSize/2)-1+jitter(2):round(stimSize/2)+eHeight+jitter(2),:) = ones([length(round(stimSize/3)-eWidth+1+jitter(1):round(stimSize/3)+eWidth+jitter(1)),length(round(stimSize/2)-1:round(stimSize/2)+eHeight),3]);
    
    % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(eMovie, 'rgb', params.meanLuminance, display);

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
