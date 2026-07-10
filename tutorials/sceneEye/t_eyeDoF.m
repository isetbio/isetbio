%% t_eyeDoF.m
%
% We suggest reading the t_eyeIntro and perhaps t_eyeFocalDistance before
% running this tutorial.
%
% This tutorial renders a retinal image of the chess set to illustrate
%
%    * Rendering with a pinhole, no lens.  Infinite depth of field
%    * Set the physiological optics and a focal distance on the rook.
%    * For a big pupil, the nearby pieces are very blurred
%    * Reduce the pupil diameter and the nearby pieces are less blurred.
%
% For speed, we do this calculation without chromatic aberration.
%
% N.B. The Chess Set scene is not part of the GitHub repository.  To
% download the PBRT scene please visit
%
%     INSTRUCTIONS WILL BE PLACED HERE.  In the mean time, ask me.
% 
% Depends on: ISET3d, ISETBio, Docker
%
% See also
%   t_eyeIntro, t_eyeFocalDistance, 

%% Initialize ISETBIO
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end
ieInit;
if ~piDockerExists, piDockerConfig; end

%% Load scene

% This version of the chess set is scaled to the right physical dimensions
thisSE = sceneEye('chessSet');

fprintf('Pupil diameter:  %0.1f mm\n',thisSE.get('pupil diameter','mm'));

%% Render a quick, low quality scene

% Through a pinhole you get a large depth of field
thisSE.set('fov',20);             % deg

thisSE.set('use pinhole',true);
thisSE.set('render type',{'radiance','depth'});
scene = thisSE.piWRS;

% Summarize what we did
thisSE.summary;

% For fun, have a look at the depth map (units are meters).  You can see
% that different chess pieces are at different distances.
scenePlot(scene,'depth map');  colorbar;

%% Big pupil case.

% Set up for optical image with lens
thisSE.set('use optics',true);

% Focus on the rook at the back (1 diopter)
thisSE.set('accommodation',1/0.6);

% Focus closer on the bishop
% thisSE.set('accommodation',1/.3);

% Upgrade the rendering quality.  
thisSE.set('pupil diameter',5);     
thisSE.set('rays per pixel',256);

% Print a summary
thisSE.summary;

% oi = thisSE.piWRS;

% Radiance only for speed
oi = thisSE.render;    
oiWindow(oi);
%{
oi = piAIdenoise(oi);
oiWindow(oi);
%}

% The overall image appearance is yellow because of the lens pigment.

%% Reducing the pupil sharpens up the nearer pieces

% Shrink the pupil.  That increases the depth of field.
thisSE.set('pupil diameter',2);   

% Shrinking the pupil also adds more rendering noise.  So we increase the
% number of rays.  This will slow the calculation.
thisSE.set('rays per pixel',512);

oi = thisSE.render('render type','radiance');    % Radiance and depth
oiWindow(oi);

thisSE.summary;

% Notice that the OI is less blurry for the chess pieces that are close. 

% There is more rendering noise, though. We reduced the pupil area by
% about a factor of 6 and only doubled the number of rays.  Impatient
% people out here in California.
%
% But we have a denoiser if you want.
%{
oi = piAIdenoise(oi);
oiWindow(oi);
%}

%% END
