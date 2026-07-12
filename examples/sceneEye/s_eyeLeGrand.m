%% s_eyeLeGrand.m
% SkipFile
% LeGrand eye examples use CPU optics renders, so keep them out of routine
% smoke-test sessions.
%
% We recommend you go through t_eyeIntro.m before running
% this tutorial.
%
% This tutorial renders the PBRT SCENE "letters at depth" using the
% LeGrand eye model.  The purpose of the script is to illustrate how to
%
%   * set up a sceneEye with the Navarro model
%   * position the camera and look at a particular scene object
%   * render with chromatic aberration
%
% Depends on: 
%    ISETBio, ISET3d, Docker
%
% Wandell, 2020
%
% See also
%   t_eye*
%

%% Check ISETBIO and initialize

ieInit;
try
    if ~piDockerExists
        piDockerConfig;
    end
catch
    warning('ISET3d not found on path. ISET3d must be your path.');
    return;
end

%% Try just the slanted edge

%{
% I need to figure out the accommodation for this case.
% It is way too blurry, as if I don't know the proper distance
%
thisSE = sceneEye('slantededge','eye model','legrand');
thisSE.set('to',[0 0 0]);

thisLight = piLightCreate('spot light 1', 'type','spot','rgb spd',[0.5 0.5 1]);
thisSE.set('light',thisLight, 'add');
thisSE.set('light',thisLight.name,'specscale',0.5);
thisSE.set('fov',12);

thisSE.set('rays per pixel',256);  % Pretty quick, but not high quality
thisSE.set('use pinhole',true);
scene = thisSE.render('docker',thisDocker);
sceneWindow(scene);

thisSE.set('use pinhole',false);
oi = thisSE.render();
oiWindow(oi);
%}

%% Here are the World positions of the letters in the scene

% The units are in meters
toA = [-0.0486     0.0100     0.5556];
toB = [  0         0.0100     0.8333];
toC = [ 0.1458     0.0100     1.6667];

%% Show the scene

% This is rendered using a pinhole so the rendering is fast.  It has
% infinite depth of field (no focal distance).
thisSE = sceneEye('letters at depth','eye model','legrand');

% Given the distance from the scene, this FOV captures everything we want
thisSE.set('fov',15);             % Degrees

% Position the eye off to the side so we can see the 3D easily
from = [0.25,0.3,-1.3];
thisSE.set('from',from);

% Look at the position with the 'B'.  The values for each of the letters
% are included above.
thisSE.set('to',toB);

% Reduce the rendering noise by using more rays. 
thisSE.set('rays per pixel',32);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('film resolution',512);  

%% Have a quick check with the pinhole
thisSE.set('use pinhole',true);

% thisSE.get('object distance')   % Default is 2.1674
% If we make it further, we can narrow the FOV, I think
% thisSE.set('object distance',6);
% thisSE.set('fov',6);

thisSE.set('render type',{'radiance','depth'});

% Summarize
thisSE.summary;

%% Render the pinhole
thisDocker = isetdocker;
scene = thisSE.render('docker',thisDocker);
scene = sceneSet(scene,'name','legrand-pinhole');
sceneWindow(scene);

% To see the depth map 
%   scenePlot(scene,'depth map');

%% Now use the optics model with chromatic aberration

% Turn off the pinhole. There is no way to adjust the LeGrand eye
% accommodation.  We can for Arizona and Navarro.
thisSE.set('use pinhole',false);

% This sets the chromaticAberrationEnabled flag and the integrator to
% spectral path.
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

% Reduce the rendering noise by using more rays. 
thisSE.set('rays per pixel',256);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('film resolution',384);     

% Summarize
thisSE.summary;

% Render
oi = thisSE.render('docker',thisDocker);
oi = oiSet(oi,'name','legrand');
oiWindow(oi);

%% END
