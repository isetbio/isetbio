%% t_slantedBarMTF.m
%
% This tutorial shows how you can render a retinal image of "slanted bar."
% We can then use this slanted bar to estimate the modulation transfer
% function of the optical system.
%
% We also show how the color fringing along the edge of the bar due to
% chromatic aberration. 
%
% We recommend you go through t_rayTracingIntroduction.m before running
% this tutorial.
%
% Depends on: pbrt2ISET, ISETBIO, Docker, ISET
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;
if ~piDockerExists, piDockerConfig; end

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
% This render takes around 2 minutes on a machine with 2 cores

% Note: The distance between the back of the lens and the front of the lens
% is approximately 7.69 mm. When we define accommodation its relative to
% the front of the lens, but for PBRT the distance to the plane is relative
% to the back of the lens. We account for this discrepancy by adding a
% slight shift. In most cases this slight difference does not make a huge
% difference, but for color fringing it does. In the future we need to fix
% this discrepancy.
planeDistance = [195 200 205] + 7.69;

for ii = 1:length(planeDistance)
    
    myScene = sceneEye('slantedBar','planeDistance',planeDistance(ii));
    myScene.name = sprintf('slantedBar_LCA_%0.2fmm',planeDistance(ii));
    
    % Zoom in to see the color fringes
    myScene.fov = 1;
    
    myScene.accommodation = 5;
    myScene.numCABands = 8;
    myScene.numRays = 64;
    myScene.resolution = 128;
    
    oi = myScene.render;
    
    vcAddObject(oi);
    oiWindow;
end


%% Calculate the MTF 
% We can use the ISO12233 standard to calculate the MTF from a slanted bar.

% First render the slanted bar. You might want to increase the numRays and
% resolution for less noisy results.
myScene = sceneEye('slantedBar','planeDistance',200);
myScene.name = 'slantedBarForMTF';
myScene.accommodation = 5;
myScene.fov = 1;
myScene.numCABands = 8;
myScene.numRays = 128;
myScene.resolution = 128;
oi = myScene.render;

% Crop the image so we only have the bar
cropRadius = myScene.resolution/(2*sqrt(2))-5;
oiCenter = myScene.resolution/2;
barOI = oiCrop(oi,round([oiCenter-cropRadius oiCenter-cropRadius ...
    cropRadius*2 cropRadius*2]));

% Convert to illuminance
% TODO: Is this the best way to go about calculating the MTF? Should we
% instead take the 550 nm photons? How should we convert from photons to
% RGB values to be passed into the ISO12233 routine? 
barImage = oiGet(barOI,'illuminance');

% Calculate MTF
deltaX_mm = oiGet(oi,'sample spacing')*10^3; % Get pixel pitch
% ISET is required for this calculation (not in ISETBIO)
% TODO: Add this in ISETBIO? 
[results, fitme, esf, h] = ISO12233(barImage, deltaX_mm(1),[1/3 1/3 1/3]);

