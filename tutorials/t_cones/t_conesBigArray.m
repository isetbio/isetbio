%% t_conesBigArray
%
%  Create a fairly large rectangular array
%
% 

%% 
ieInit

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

% Using the pull down in the window, you can choose the image that is
% displayed. You can bring up one or the other images with the 'show'
% argument
cMosaic.window('show','mean absorptions');

%% Or the cone mosaic
cMosaic.window('show','cone mosaic');

%% No photocurrent because we didn't compute it

% But if we add eye movements, recompute, then we can show
%{
   cMosaic.emGenSequence(100);
   cMosaic.compute(oi);
   cMosaic.computeCurrent;
   cMosaic.window('show','mean photocurrent');
%}
%% END