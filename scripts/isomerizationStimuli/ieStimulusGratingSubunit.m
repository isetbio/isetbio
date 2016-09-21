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
addParameter(p,'freq',          5,  @isnumeric);
addParameter(p,'radius',          0,  @isnumeric);
addParameter(p,'theta',          0,  @isnumeric);
addParameter(p,'side',          0,  @ischar);
addParameter(p,'downsamplefactor', 1,  @isnumeric);
p.parse(varargin{:});
params = p.Results;
fov = params.fov;

wFlag = ieSessionGet('wait bar');

%% Compute a Gabor patch scene as a placeholder for the bar image

% Create display
display = displayCreate;

% Set up scene, oi and sensor
scene = sceneCreate();
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

%% Initialize the optics and the sensor
oi  = oiCreate('wvf human');

if params.radius == 0
    % Assume foveal cone density
    absorptions = sensorCreate('human');

else
    % Generate cone mosaic with density according to eccentricitiy
    coneP = coneCreate; % The default cone properties
    retinalRadiusDegrees = params.radius;
    retinalPolarDegrees  = params.theta;
    whichEye             = params.side;
    retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; 
    absorptions = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
end

% Get the cone mosaic and plot it
cone_mosaic = sensorGet(absorptions,'cone type');
% vcNewGraphWin; imagesc(cone_mosaic); colormap jet;

% Set cone mosaic
cone_mosaic = 3*ones(sensorGet(absorptions,'size'));
% cone_mosaic = 2+round(rand(sensorGet(absorptions,'size')));
% blueIndices = find(cone_mosaic==0);
% cone_mosaic(blueIndices) = 2;
absorptions = sensorSet(absorptions,'cone type',cone_mosaic);

% Scale the size of cone mosaic approrpiately
absorptions = sensorSetSizeToFOV(absorptions, params.fov, scene, oi);

% Set aspect ratio
% At this point, we have the right cone density and the right number in
% cols, now we just need to set the rows to have the same aspect ratio as
% the input movie.
sceneSize = sceneGet(scene,'size');
sensorSize = sensorGet(absorptions,'size');
aspectRatioMovie = sceneSize(1)/sceneSize(2);
absorptions = sensorSet(absorptions,'size',[aspectRatioMovie*sensorSize(2) sensorSize(2)]);

absorptions = sensorSet(absorptions, 'exp time', params.expTime); 
absorptions = sensorSet(absorptions, 'time interval', params.timeInterval); 

%% Compute a dynamic set of cone absorptions for moving bar

fprintf('Computing cone isomerizations:    \n');

% ieSessionSet('wait bar',true);
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
% The number of steps must be smaller than the width of the scene
nSteps = params.nSteps;

downSampleFactor = 8;

for t = 1 : nSteps
    if wFlag, waitbar(t/nSteps,wbar); end
        
%     sceneSize = sceneGet(scene, 'size');
    sceneSize = [params.row params.col];
    barMovie = ones([sceneSize(1)+2*params.barWidth,sceneSize(2), 3])*0.001;  % Gray background
    
    %horz
%     nStripes = floor((sceneSize(1)+params.barWidth)/params.barWidth);
    % vert
    nStripes = floor((sceneSize(2)+params.barWidth)/params.barWidth);
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
    randWalk(t)=floor(t/2);
    for stripeInd = 1:2:nStripes+2
        % barMovie(randWalk(t)+(stripeInd-1)*params.barWidth+params.barWidth:randWalk(t)+(stripeInd)*params.barWidth+params.barWidth-1,:,:) = 1;          % White bar
        
        % horizontal stripes
%         barMovie((stripeInd-1)*params.barWidth+params.barWidth:(stripeInd)*params.barWidth+params.barWidth-1,:,:) = .75;          % White bar
        
        % vertical stripes
        barMovie(:,(stripeInd-1)*params.barWidth+params.barWidth:(stripeInd)*params.barWidth+params.barWidth-1,:) = 1;          % White bar
    end
%     barMovie = circshift(barMovie,randWalk(t),2);
%     if mod(floor(t/20),2) == 0
%         barMovie = 1-barMovie;
%     end
    
%     barMovie = 0.5+1*(barMovie - 0.5)*sin(2*pi*((t-1)/200));

    % barMovie = 0.5+1*(barMovie - 0.5)*sin(2*pi*((t-1)/(100/params.freq)));
    barMovie = ((0.5+1*(barMovie - 0.5)*(sin(2*pi*((t-1)/(100/params.freq)))))>0.5);
      
%     barMovieResize = barMovie(params.barWidth+1:params.barWidth+sceneSize(1),1:sceneSize(2),:);
    barMovieResize = barMovie(params.barWidth+1:params.barWidth+sceneSize(1),1:sceneSize(2),:);
%     if t < 6
%         barMovieResize = ones(size(barMovieResize))*0.5;
%         barMovieResize(1,1,:) = ones(1,1,3);
%     end
%     barMovieResize(1,1,:) = ones(1,1,3);
%     barMovieResize(1,2,:) = zeros(1,1,3);
    
    barMovie(1,1,:) = ones(1,1,3);
    barMovie(1,2,:) = zeros(1,1,3);
    % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);

    scene = sceneSet(scene, 'h fov', fov);
    if t ==1
        % sceneRGB = zeros([sceneGet(scene, 'size'), nSteps, 3]);
        sceneRGB = zeros([size(barMovie,1) size(barMovie,2) nSteps, 3]);
    end
    
    % Get scene RGB data    
%     sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
%     sceneRGB(1,1,t,:) = sceneRGB(2,1,t,:);
%     sceneRGB(1,2,t,:) = sceneRGB(2,2,t,:);
    
    sceneRGB(:,:,t,:) = barMovie;
    
% % % % % % % % % % % % % % % % % % % % %     
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    for tNew = 1:downSampleFactor
    
        % Compute absorptions
        
        cone_mosaic = 3*ones(sensorGet(absorptions,'size'));
        % cone_mosaic = 2+round(rand(sensorGet(absorptions,'size')));
        % blueIndices = find(cone_mosaic==0);
        % cone_mosaic(blueIndices) = 2;
        absorptions = sensorSet(absorptions,'cone type',cone_mosaic);
        absorptions = sensorComputeNoiseFree(absorptions, oi);
        
        if t == 1
            volts = zeros([sensorGet(absorptions, 'size') params.nSteps]);
        end
        
        volts(:,:,downSampleFactor*(t-1)+tNew) = sensorGet(absorptions, 'volts');
    end
% % % % % % % % % % % % % % % % % % % % % % % % % 


    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end


% % Compute optical image
% oi = oiCompute(oi, scene);
% 
% % Compute absorptions
% absorptions = sensorCompute(absorptions, oi);
% 
% volts(:,:,1) = sensorGet(absorptions, 'volts');

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
