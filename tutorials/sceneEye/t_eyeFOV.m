%% t_eyeFOV
%
% Field of view calculations with the sceneEye model.  This illustrates the
% relationship between the field of view and the 'semi diam'.
%
% The 'semi diam' is the semi diameter of the curved retinal surface.  See
% one of the PPT talks (e.g. Spitschan OSA seminar) for definition graph.
%
% See also
%  t_eye*
%

%% Check ISETBIO and initialize

% The sceneEye modeling uses ISETBio and ISET3d.  So we check.
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end

% Then we initialize ISETBio and make sure the user has Docker configured.
ieInit;
if ~piDockerExists, piDockerConfig; end

%% Load a scene and set an eye model

% thisEye = sceneEye('slantedbar');
% thisEye.summary;

% There are a number of PBRT scenes stored in piRootPath/data/V3.  This is
% one of them.  To read in the recipe for this scene, we use this command.
thisEye = sceneEye('chessset');

% thisEye.set('mmUnits',true);
thisEye.summary;

% thisR = piRecipeDefault('scene name','sanmiguel');

%%  The original scene

% This has a large field of view.

oi = thisEye.render;
oiWindow(oi);

%% Narrow the FOV and move the camera
thisEye.set('fov',2);
thisEye.set('rays per pixel',128);

% Move the camera in the negative direction (away from the scene)
thisFrom = thisEye.get('from');
thisEye.set('from',thisFrom + [0 0 -3]);
oDist = thisEye.get('object distance','m');

thisEye.set('focal distance',oDist);
oi = thisEye.render;
oiWindow(oi);

%% THe relationship between the field of view and semidiam

thisEye.get('fov')
thisEye.get('retina semidiam','mm')

thisEye.set('fov',10);
thisEye.get('retina semidiam','mm')

%% Set the fov 
thisEye.set('fov',2);
thisEye.get('retina semidiam','mm')

thisEye.get('fov')

%% Now set the semidiam
thisEye.set('retina semidiam',4);
thisEye.get('fov')

thisEye.set('retina semidiam',2);
thisEye.get('fov')

thisEye.set('retina semidiam',1);
thisEye.get('fov')

%%
thisEye.set('chromatic aberration',true)
oi = thisEye.render;
oiWindow(oi);

thisEye.set('chromatic aberration',false)
thisEye.summary
oi = thisEye.render;
oiWindow(oi);

%% END