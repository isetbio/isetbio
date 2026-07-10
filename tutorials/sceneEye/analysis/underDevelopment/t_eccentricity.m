% SkipFile
%% t_eccentricity.m
%
% This script demonstrates the eccentricity parameter in the sceneEye
% object of ISETBIO. 
% 
% When eccentricity is set, the center of the rendered image will
% correspond to the given x/y eccentricity on the retina. 
%
% "Developer" note: 
% The way this is done in sceneEye.render is a bit of a hack, in that we
% will instruct pbrt-v3-spectral to render a very large retinal image that
% encompasses the desired eccentricity, but then also set a crop window
% which will only render the section of the image we are interested in.
% However, these values are all calculated automatically and for the most
% part this isn't something that a regular user needs to worry about. 
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;

%% Render the original scene
% By default, the center of the image corresponds to the center of the
% retina.

myScene = sceneEye('numbersAtDepth');
myScene.fov = 30;

myScene.name = '0_0_ecc';

myScene.numRays = 64;
myScene.resolution = 106;
oi = myScene.render;

ieAddObject(oi);
oiWindow;

%% Render a couple of smaller images a different locations on the retina

% TODO: There seems to be some issues with the magnification with the
% eccentricity calculations. We get a slight decrease in magnification if
% we render with a larger FOV and then crop to the smaller size vs
% rendering directly with the smaller size.
% For example, say we have two renders...
%
% Render 1: 30 degrees FOV with 64 pixel width (centered at 0,0)
% Render 2: 47.91 degrees FOV with 106 pixel width (centered at 0,0)
%           Cropped to center 64x64 pixels
%
% Render 1 and 2 will be off in magnfication by a few pixels.
%
% I believe all the math is correct, it is still unknown why this is
% happening. Does it have to do with the spherical nature of the retina and
% how we sample it? I think this might be the case.
% Needs some deeper investigation.

myScene = sceneEye('numbersAtDepth');
myScene.fov = 30;

ecc = [10 0;
    10 10;
    -10 -10];

myScene.resolution = 64;

for ii = 1:length(ecc)
    
    currEcc = ecc(ii,:);
    myScene.eccentricity = currEcc;
    myScene.name = sprintf('%i_%i_ecc',currEcc(1),currEcc(2));
    
    oi = myScene.render;
    ieAddObject(oi);
    oiWindow;
    
end
