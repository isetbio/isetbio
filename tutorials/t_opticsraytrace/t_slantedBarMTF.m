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
    
% TODO:
% Why is the slanted edge off center? Is the plane not centered correctly?

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

%% Turn on chromatic aberration to show color fringing.
% We can render chromatic aberration in the eye by tracing one ray per band
% of wavelength. The parameter, numCABands determines the number of band we
% will sample. We will trace a total of numRay x numCABands rays, meaning
% that the rendering will be ~(numCABands) times slower.
% 
% As we move the plane in and out of focus we can see the color fringes
% change due to longitudinal chromatic aberration.
%
% This render takes around 4-5 minutes on a machine with 2 cores

% Note: The distance between the back of the retina and the front of the
% lens is approximately 7.69 mm. When we define accommodation its relative
% to the front of the lens, but for PBRT the distance to the plane is
% relative to the retina. We account for this discrepancy by adding a
% slight shift. In most cases this slight difference does not make a huge
% difference, but for color fringing it does. In the future we need to fix
% this discrepancy. 
planeDistance = [95 100 110] + 7.69;

for ii = 1:length(planeDistance)
    
    myScene = sceneEye('slantedBar','planeDistance',planeDistance(ii));
    myScene.name = sprintf('slantedBar_LCA_%0.2fmm',planeDistance(ii));
    
    myScene.accommodation = 10;
    myScene.numCABands = 8;
    myScene.numRays = 128;
    myScene.resolution = 128;
    
    oi = myScene.render;
    
    vcAddObject(oi);
    oiWindow;
end

% TODO:
% Show color fringing from longitudinal chromatic aberration
% Show how you can use ISO12233 to estimate the MTF

%% Calculate the MTF 
% We can use the ISO12233 standard to

