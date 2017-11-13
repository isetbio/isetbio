%% t_slantedBarMTF.m
%
% This tutorial shows how you can render a slanted bar through the eye. We
% can then use this slanted bar to estimate the modulation transfer
% function of the optical system. 
%
% We also show how the MTF changes with accommodation, as well as the color
% fringing from the longitudinal chromatic aberration. 
%
% We recommend you go through t_rayTracingIntroduction.m before running
% this tutorial.
%
% Depends on: pbrt2ISET, ISETBIO, Docker
%
% TL ISETBIO Team, 2017
    

%% Initialize ISETBIO
ieInit;

%% Render a fast image of the slanted bar first

% The slanted bar scene consists of a square plane (1x1 m) that is
% split in half diagonally. The bottom left half is white while the top
% right half is black. By default the plane is placed at [0 1 0] meters,
% but we can change that by given sceneEye an optional 'planeDistance'
% input. 
myScene = sceneEye('slantedBar','planeDistance',500); % Create a slanted bar at 0.5 meter
myScene.name = 'slantedBarFast';
myScene.numRays = 64;
myScene.resolution = 128; 

oi = myScene.render;

vcAddObject(oi);
oiWindow;

%% Try moving the slanted bar in and out of focus

planeDistance = [100 300 500 1000];

for ii = 1:length(planeDistance)
    
    myScene = sceneEye('slantedBar','planeDistance',planeDistance(ii));
    myScene.name = sprintf('slantedBar_%0.2fmm',planeDistance(ii));
    myScene.numRays = 64;
    myScene.resolution = 128;
    myScene.accommodation = 2; 
    
    oi = myScene.render;
    
    vcAddObject(oi);
    oiWindow;

end

%%
% TODO:
% Show color fringing from longitudinal chromatic aberration
% Show how you can use ISO12233 to estimate the MTF



