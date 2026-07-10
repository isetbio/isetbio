%% t_eyeNavarro.m
%
% We recommend you go through t_eyeIntro.m before running
% this tutorial.
%
% This tutorial renders the PBRT SCENE "letters at depth" using the Navarro
% eye model.  The script illustrates how to
%
%   * set up a sceneEye with the Navarro model
%   * position the camera to center on a specific scene object
%   * update the Navarro lens files for different accommodations
%
% Depends on: 
%    ISETBio, ISET3d, Docker
%
% Wandell, 2020
%
% See also
%   t_eyeArizona, t_eyeLeGrand
%

%% Check ISETBIO and initialize

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
thisSE = sceneEye('letters at depth');

% First render with pinhole
thisSE.set('use pinhole',true);

thisSE.set('render type',{'radiance','depth'});

% Position the eye off to the side so we can see the 3D easily
from = [0.25,0.3,-0.2];
thisSE.set('from',from);

% Look at the position with the 'B'.
thisSE.set('to',toB);

% Set its distance
thisSE.set('object distance',1);  % meters

thisSE.set('to',toA); distA = thisSE.get('object distance');
thisSE.set('to',toB); distB = thisSE.get('object distance');
thisSE.set('to',toC); distC = thisSE.get('object distance');
thisSE.set('to',toB);

% To just see the 'B' at higher resolution use a small FOV
thisSE.set('fov',40);

% Render the scene
thisSE.set('render type', {'radiance','depth'});

%% Render as a scene with isetdocker

thisDocker = isetdocker;
scene = thisSE.piWRS('docker',thisDocker,'name','pinhole');

% You can see the depth map if you like
%   scenePlot(scene,'depth map');

%% Now use the optics model with chromatic aberration

% Use the model eye
thisSE.set('use optics',true);

thisSE.set('pupil diameter',3);

% True by default anyway
% thisSE.set('mmUnits', false);

% We turn on chromatic aberration.  That changes the lens and IOR files
% that piWrite produces for the optical model.

% This sets the chromaticAberrationEnabled flag and the integrator to
% spectral path.
% Now works in V4 - May 28, 2023 (ZL)
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

% We can reduce the rendering noise by using more rays. This takes a while.
thisSE.set('rays per pixel',256);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('spatial samples',384);     

% Ray bounces
thisSE.set('n bounces',3);

% We want the scene to be around 5-10 deg so we do not need a lot of
% samples to resolve the blur.  This renders only the region around
% the 'B'.
thisSE.set('fov',7);             % Degrees

sampleSpacing = thisSE.get('sample spacing');
if ~isempty(sampleSpacing)
    fprintf('Sample spacing: %.3f um\n',sampleSpacing*1e3);
end
%% Change the accommodation.  But look at 'B'.

% Focus on the A
thisSE.set('accommodation',1/distA);  

% Summarize
thisSE.summary;


%% Write the Navarro lens model for this accommodation

% The pinhole render above is reliable and fast.  The optical sceneEye
% rendering path still needs work in isetdocker, so here we write the
% updated Navarro lens and IOR files without rendering.
piWrite(thisSE.recipe);
fprintf('Accommodation A lens file: %s\n',thisSE.recipe.get('lens file'));

%% Change the accommodation.  But look at 'B'.

% Focus on the B
thisSE.set('accommodation',1/distB);  

% Summarize
thisSE.summary;

piWrite(thisSE.recipe);
fprintf('Accommodation B lens file: %s\n',thisSE.recipe.get('lens file'));

%% Set accommodation to a different distance.

% Focus on the C
thisSE.set('accommodation',1/distC);  

thisSE.summary;
piWrite(thisSE.recipe);
fprintf('Accommodation C lens file: %s\n',thisSE.recipe.get('lens file'));

%% END
