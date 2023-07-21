%% Compute with iset3d and a cMosaic from ISETBio
%
% See also
%   ISET3d and sceneEye
%

%%
ieInit;
if ~piDockerExists, piDockerConfig; end

%% Read and set the recipe
thisSE = sceneEye('chessSet','eye model','navarro');
thisSE.set('object distance',1);
thisSE.set('use pinhole',true);
thisSE.set('rays per pixel',128);  % Pretty quick, but not high quality
thisSE.summary;
thisSE.piWRS; 

%% Cone mosaic

thisSE.set('use pinhole',false);
thisSE.set('object distance',5);
thisSE.set('fov',10);

dockerWrapper.reset;
thisD = dockerWrapper.humanEyeDocker;
oi = thisSE.piWRS('docker wrapper',thisD);

%% Render on a cone mosaic

% This is a lot of cones and it takes a while (2 min?)
cm = cMosaic('sizeDegs',[6,3],'eccentricityDegs',[2,0]);

% The high density around the fovea is obvious.
cm.plot('mosaic');

%% Compute 

% The spatial resolution of the rendered image is pretty low, so you get a
% warning message.
excitations = cm.compute(oi);

cm.plot('excitations',excitations);

%% END