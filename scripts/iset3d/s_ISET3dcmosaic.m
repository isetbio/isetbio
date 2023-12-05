%% Integrating iset3d and isetbio
%
% These calculations show how to set up a scene and calculate the cone
% excitations for a 3D scene and a human eye model.
%
% We load the Chess Set Pieces and a human eye model (navarro).
%
% We can adjust the viewing situation, the size of the retina, and then
% calculate a cone mosaic and cone excitations
%
% See also
%   ISET3d/tutorials/sceneEye/t_eye*

%%
ieInit;
if ~piDockerExists, piDockerConfig; end

%% Read the recipe
thisSE = sceneEye('chess set','eye model','navarro');

%% Rendering parameters

thisSE.set('rays per pixel',32);  % 32 is Pretty quick, but not high quality

thisSE.set('fov',40);

% Render the scene
thisSE.set('render type', {'radiance','depth'});

thisSE.set('use pinhole',true);

thisDocker = dockerWrapper;
scene = thisSE.piWRS('docker wrapper',thisDocker,'name','pinhole');

%% Eye parameters:  position, focal distance and field of view

thisSE.set('use optics',true);

thisSE.set('pupil diameter',3);

% thisSE.set('to',[-0.09 0.05 0.01]);  % Puts the small pawn in the center
thisSE.set('to',[-0.04 0.05 0.01]);    % Puts the small pawn in the center
thisSE.set('object distance',1);       % Distance between 'from' to 'to'
thisSE.set('focal distance',0.8);      % Focus just in front if the 'to'
% piAssetGeometry(thisSE.recipe);

thisSE.set('retina Semi Diam',1)     % This controls the size of the retina

thisSE.set('spatial samples',[320 320]*2);  % 

%%
thisSE.lensDensity = 1;
thisDocker = dockerWrapper.humanEyeDocker;
thisSE.piWRS('docker wrapper',thisDocker,'name','Density 1');

%%
thisSE.lensDensity = 0;
thisSE.piWRS('docker wrapper',thisDocker,'name','Density 0');

%% Cone mosaic

% Pick a region of the image.  We need to add some select options to the oi
% window (and scene window).
cm = cMosaic('sizeDegs',[6,4],'eccentricityDegs',[-3,0]);
cm.visualize();

%% Compute 

cm.integrationTime = 0.05;
excitations = cm.compute(oiLeft);

params = cm.visualize('params');

params.activation = excitations.^0.5;
params.activationColorMap = hot(1024);
cm.visualize(params);

%%  Plot the histogram of excitations

r = excitations(cm.lConeIndices);
g = excitations(cm.mConeIndices);
s = excitations(cm.sConeIndices);

%%  Histogram of LMS excitations

ieNewGraphWin;
histogram(r(:),'FaceColor','r'); hold on;
histogram(g(:),'FaceColor','g'); hold on;
histogram(s(:),'FaceColor','b'); hold off;
grid on;
xlabel('Excitations'); ylabel('N cones');
legend({'L-cones','M-cones','S-cones'})
