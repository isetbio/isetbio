%% t_eyeArizona
% SkipFile
% Uses obsolete dockerWrapper setup; needs sceneEye rewrite.
%
% Run the Arizona eye model
%
% The script runs the Arizona eye model to show it runs.  But it contains
% comments to show how to turn on chromatic aberration, narrow the FOV, and
% look at the spread in more detail.
%
% See also
%   t_eyeNavarro, t_eyeLeGrand

%% Initialize paths and such

ieInit;
if ~piDockerExists, piDockerConfig; end

%% Here are the World positions of the letters in the scene

% The units are in meters
toA = [-0.0486     0.0100     0.5556];
toB = [  0         0.0100     0.8333];
toC = [ 0.1458     0.0100     1.6667];

%% Show the scene

% This is rendered using a pinhole so the rendering is fast.  It has
% infinite depth of field (no focal distance).
thisSE = sceneEye('letters at depth','eye model','arizona');

% Position the eye off to the side so we can see the 3D easily
from = [0.25,0.3,-0.2];
thisSE.set('from',from);

% Look at the position with the 'B'.  The values for each of the letters
% are included above.
thisSE.set('to',toB);

% Have a quick check with the pinhole
thisSE.set('use pinhole',true);

% Given the distance from the scene, this FOV captures everything we want
thisSE.set('fov',30);             % Degrees

%%  Render

% Summary of status
thisSE.summary;

thisDocker = isetdocker;
thisSE.piWRS('docker',thisDocker,'name','pinhole');

%% Now use the optics model with chromatic aberration

% Turn off the pinhole.  The model eye (by default) is the Navarro model.
thisSE.set('use optics',true);

% True by default anyway
thisSE.set('mmUnits', false);

% We turn on chromatic aberration.  That slows down the calculation, but
% makes it more accurate and interesting.  We often use only 8 spectral
% bands for speed and to get a rought sense. You can use up to 31.  It is
% slow, but that's what we do here because we are only rendering once. When
% the GPU work is completed, this will be fast!

% This sets the chromaticAberrationEnabled flag and the integrator to
% spectral path.
% Now works in V4 - May 28, 2023 (ZL)
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

% Distance in meters to objects to govern accommodation.
thisSE.set('to',toA); distA = thisSE.get('object distance');
thisSE.set('to',toB); distB = thisSE.get('object distance');
thisSE.set('to',toC); distC = thisSE.get('object distance');
thisSE.set('to',toB);

% We can reduce the rendering noise by using more rays. Sometimes we
% use 512 instead of 256.
thisSE.set('rays per pixel',256);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('spatial samples',256);     

% Ray bounces
thisSE.set('n bounces',3);

%% Accommodate to letter A distance (in diopters)

thisSE.set('accommodation',1/distA);
% thisSE.get('accommodation')

% Summarize
thisSE.summary;

%%
thisSE.piWRS('docker',thisDocker,'name','arizona-A');

%% Make an oi of the chess set scene using the LeGrand eye model

thisSE.set('accommodation',1/distC);  
thisSE.piWRS('docker',thisDocker,'name','arizona-C');

%% Have a look with the slanted bar scene

% Commented out because it takes a while to run.  But in a way, seeing the
% chromatic aberration is the point.  So, I put it in here.  The slanted
% bar is at the focal distance.
thisSE = sceneEye('slantedEdge','eye model','arizona');
thisSE.set('to',[0 0 0]);

thisLight = piLightCreate('spot light 1', 'type','spot','rgb spd',[0.5 0.5 1]);
thisSE.set('light',thisLight, 'add');
thisSE.set('light',thisLight.name,'specscale',0.5);
thisSE.set('fov',2);

thisSE.set('rays per pixel',256);  % Pretty quick, but not high quality

thisSE.set('use pinhole',true);
thisSE.piWRS('docker',thisDocker);  % Render and show

%% Now the human eye model

thisSE.set('use pinhole',false);

% This sets the chromaticAberrationEnabled flag and the integrator to
% spectral path.
% Now works in V4 - May 28, 2023 (ZL)
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

thisSE.set('object distance',20);

oi = thisSE.piWRS('docker',thisDocker,'show',false,'name','8bands');

% Maybe we should set a denoise flag in piWRS?
oi = piAIdenoise(oi); 

oiWindow(oi);


thisSE.summary;

%% END
