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
s = sceneSet(s, 'fov', 1);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

%% Build a default cone mosaic and compute isomerizatoins

% Get the default set of cone mosaic parameters.
[~,cmParams] = cMosaic('params');

% Generate off-axis mosaic
cmParams.eccentricityDegs = [0 0];   % The size of the cones out here is about the same as the mouse?
cmParams.sizeDegs = [1.5 1.5];
cmParams.micronsPerDegree = oiGet(oi,'distance per degree','um');
cm = cMosaic(cmParams);

cm.visualize;

%% Compute isomerizations for each eye position.
[noiseFree, noisy] = cm.compute(oi);

vParams = cm.visualize('params');
vParams.activation = noisy;
vParams.activationColorMap = gray(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);

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
% cMosaic.plot('meanabsorptions');

%% END
