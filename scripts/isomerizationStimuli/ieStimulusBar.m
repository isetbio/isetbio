function iStim = ieStimulusBar(varargin)
%% Compute a dynamic set of cone absorptions for moving bar
%
%
% We want to produce a scene video that translates into an oi video that
% becomes a cone absorption video.  At present coneAbsorptions ONLY does
% this using eye movements, not by creating a series of images.  This code
% represents our first effort to produce dynamic scenes.
%
% We are literally going to recreate a set of scenes with different phase
% positions and produce the scenes, ois, and cone absorptions by the loop.
% The result will be a time series of the cone photon absorptions.
%
% We are reluctant to make scene(:,:,:,t) because we are frightened about
% the size.  But it still might be the right thing to do.  So the code here
% is an experiment and we aren't sure how it will go.

% sceneRGB = zeros([sceneGet(scene, 'size') params.nSteps 3]); % 3 is for R, G, B
% sensorPhotons = zeros([sensorGet(sensor, 'size') params.nSteps]);
% stimulus = zeros(1, params.nSteps);
fprintf('Computing cone isomerization:    \n');

% ieSessionSet('wait bar',true);
wFlag = ieSessionGet('wait bar');
if wFlag, wbar = waitbar(0,'Stimulus movie'); end

% Loop through frames to build movie
for t = 1 : params.nSteps
    if wFlag, waitbar(t/params.nSteps,wbar); end
        

%     if t == 1 
        barMovie = zeros(params.row,params.col,3);
%     end
    barMovie(:,1+t:t+params.barwidth,:) = 0.5 + 0.499*ones(params.row,params.barwidth,3);
    barMovie(1,1,:) = 0;

    % % % % Generate scene object from stimulus RGB matrix and display object
    scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);

    scene = sceneSet(scene, 'h fov', fov);

    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);    
    
    % Compute absorptions
    sensor = sensorCompute(sensor, oi);

    if t == 1
        volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
    end
    
    volts(:,:,t) = sensorGet(sensor, 'volts');
    
    % vcAddObject(scene); sceneWindow
end

if wFlag, delete(wbar); end

% Set the stimuls into the sensor object
sensor = sensorSet(sensor, 'volts', volts);
% vcAddObject(sensor); sensorWindow;

% These are both the results and the objects needed to recreate this
% script. So calling isomerizationBar(iStim) should produce the same
% results
%
iStim.oi      = oi;
iStim.display = display;
iStim.params  = params;
iStim.sensor  = sensor;
end
