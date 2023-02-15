%% Introduction to the cone mosaic object.
%
% Description:
%    Create a rectangular cone mosaic object and compute cone
%    isomerizatoins across a set of small eye movements.
%
%    Visualize the results in coneMosaic window.
%
% See also
%   t_cones*

%% Initialize and clear
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 2.5);
sceneWindow(s);

%% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);
oiWindow(oi);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show about half the scene. Speeds things up.
cMosaic.setSizeToFOV(0.9 * sceneGet(s, 'fov'));

% You can see the field of view for this cone mosaic object, along with
% other parameters, within the coneMosaic object:
% cMosaic.fov;

%% Generate a sequence of eye posistions.

% This function creates an eye movement object (see t_fixationalEM.m) and
% automatically generates a path for this given cone mosaic.
cMosaic.emGenSequence(50);

%% Create a shadow pattern 

sz = oiGet(oi,'size');
shadow = imread(fullfile(isetbioRootPath,'local','bloodVessels.png'));
shadow = sum(shadow,3);
shadow = ieScale(shadow,0,.4);
shadow = 1 - shadow;
g = fspecial('gaussian',[9,9]);
shadow = convolvecirc(shadow,g);

%% Compute isomerizations for each eye position.

% We want the shadow to be attached to the cones (not the OI).  It
% should be inserted here, as part of the compute
cMosaic.compute(oi);

shadow = imresize(shadow,[cMosaic.rows,cMosaic.cols]);
cMosaic.absorptions = cMosaic.absorptions .* shadow;

%% Bring up a window so that we can look at things.
%
% Using the pull down in the window, you can look at
% the mosaic, the isomerizations for one fixation, or
% the movie of fixations.
cMosaic.window;

%% Plot things individually
% Instead of using the cone mosaic window, you can call the plot function
% directly. 

% cMosaic.plot('Cone mosaic');
% cMosaic.plot('mean absorptions');

%% END
