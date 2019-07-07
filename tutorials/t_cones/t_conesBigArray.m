%% t_conesBigArray
%
%  Create a fairly large rectangular array
%
% 

%% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 3);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show about half the scene. Speeds things up.
cMosaic.setSizeToFOV(sceneGet(s, 'fov'));

%% Compute isomerizations for each eye position.
cMosaic.compute(oi);

%% Bring up a window so that we can look at things.
%
% Using the pull down in the window, you can look at
% the mosaic, the isomerizations for one fixation, or
% the movie of fixations.
cMosaic.window;

%%