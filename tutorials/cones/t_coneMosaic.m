%% t_coneMosaic
%
% Introduction to the cone mosaic object.

%% Initialize and clear
ieInit

%% Build a scene and oi for computing
%
% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s,'fov',1);
ieAddObject(s);

%% Then the oi
oi = oiCreate;
oi = oiCompute(oi,s);
ieAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute from the OI
cMosaic = coneMosaic;  % Create the object
cMosaic.setSizeToFOV(0.5*sceneGet(s,'fov'));   % Show about half the scene
cMosaic.emGenSequence(100);
cMosaic.compute(oi);
cMosaic.window;
