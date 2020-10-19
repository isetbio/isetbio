%% Show a letter on a display imaged at the level of the cone mosaic
%
%
% TODO:
%   Write a routine that takes the outline of the letter and superimposes
%   it on the color image of the cones.
%
%   Implement this for the hex mosaic.
%
% Wandell, 2019
%
% See also
%   sceneCreate('letter', ...), displayCreate, ....
%

%%
ieInit;

%% Create a letter on a display

% family, size, dpi
font = fontCreate('A', 'Georgia', 14, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',0.6);

% We should pad the scene so the eye movements do not move the scene beyond
% the array

% Here is the scene
sceneWindow(scene);

%% Push the scene through human optics

oi = oiCreate;
oi = oiCompute(oi,scene);
oiWindow(oi);
%%  Now image it on the cone mosaic with some fixational eye movements

cones = coneMosaic;
cones.setSizeToFOV(1.3*sceneGet(scene,'fov'));
cones.emGenSequence(50);
cones.compute(oi);
cones.window;

%% For a hex mosaic now

% Set up the parameters and make this version work, next.
%
%{
 resampleFactor = 4;
 conesH = coneMosaicHex(resampleFactor,'fovDegs',0.5);
 conesH.emGenSequence(50);
 conesH.compute(oi);
 conesH.window;
%}
