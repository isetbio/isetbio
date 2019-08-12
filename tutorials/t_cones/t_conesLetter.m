%% Show a letter on a display imaged at the level of the cone mosaic
%
%

%%
ieInit;

%% Create a letter on a display

% family, size, dpi
font = fontCreate('A', 'Georgia', 24, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',0.5);

% Here is the scene
sceneWindow(scene);

%% Push the scene through human optics

oi = oiCreate;
oi = oiCompute(oi,scene);

%%  Now image it on the cone mosaic with some fixational eye movements

cones = coneMosaic;
cones.setSizeToFOV(sceneGet(scene,'fov')*1.5);
cones.emGenSequence(100);
cones.compute(oi);
cones.window;

%% For a hex mosaic now

%%
