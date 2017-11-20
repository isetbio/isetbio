%% t_planarImage.m
%
% This tutorial demonstrates how to use the textured plane scene. With this
% scene, we can load up an image texture, place it on a scaled, 2D plane at
% a certain distance away from the eye. We can then render the image
% through the eye to obtain the retinal image.
% 
% This can be used as an approximation of the eye viewing a flat "display."
% It's not quite the same as using a simulated display in ISETBIO with real
% primaries, but it does allow us to use the 3D eye model with a flat
% image.
%
% For example, you can use this textured plane to approximate the following steps: 
% 1. Render a retinal image from a 3D scene,
% 2. Apply some processing on the resulting retinal image
% 3. Place it on a flat display of size xxx at yyy distance away from the
% viewer.
% 4. Render a new retinal image.
%
% Depends on: pbrt2ISET, ISETBIO, Docker, ISET
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;

%% Load up the textured plane scene

% We're going to place the "display" 1 meter away from the eye with a 40
% degree FOV.
distance = 1000;
planeFOV = 40;
width = 2*tand(planeFOV/2)*distance;
size = [width width];

% We're going to use the following texture on the plane
imageTexture = fullfile(piRootPath,'data','imageTextures','squareResolutionChart.exr');

myScene = sceneEye('texturedPlane',...
                   'planeDistance',distance,...
                   'planeSize',size,...
                   'planeTexture',imageTexture);
               
%% Change eye properties

% We can have the display image approximately cover the entire retinal
% image by matching the FOV. 
myScene.fov = planeFOV;

% Accommodate to the plane
myScene.accommodation = 1/(distance*10^-3); 


%% Render
oi = myScene.render;
vcAddAndSelectObject(oi);
oiWindow


                
                
