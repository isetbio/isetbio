%% Introduction to the cone mosaic object.
%
% Description:
%    Show how to create a cone mosaic object and compute cone
%    isomerizatoins across a set of small eye movements. This lets you look
%    at the result in the coneMosaic window.
%
% Notes:
%    * [Note: DHB - Either we should show how to access fields of the
%      mosaic programatically here, or we should point to a different
%      tutorial that does so.]
%

%% Initialize and clear
clear;
close all;
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 1);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show about half the scene. Speeds things up.
cMosaic.setSizeToFOV(0.1 * sceneGet(s, 'fov'));

%% Generate a sequence of 100 eye posistions.
cMosaic.emGenSequence(100);

%% Compute isomerizations for each eye position.
cMosaic.compute(oi);

%% Bring up a window so that we can look at things.
%
% Using the pull down in the window, you can look at
% the mosaic, the isomerizations for one fixation, or
% the movie of fixations.
cMosaic.window;

%%
