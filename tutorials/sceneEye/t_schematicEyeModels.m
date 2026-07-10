%% t_schematicEyeModels
% SkipFile
% Multi-model sceneEye comparison with several CPU human-eye renders. Keep it
% interactive rather than running it in the tutorial smoke test.
% Render the same scene using a couple of different eye models.

%% Initialize
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end
ieInit;
clear; close all;

%% Load up a scene

thisSE = sceneEye('numbersAtDepth','human eye','Navarro');
thisSE.set('mmUnits', true);

%{
% Three letters on a checkerboard background. A is at 1.4 dpt, B is at 1
% dpt, and C is at 0.6 dpt. 
thisSE = sceneEye('lettersAtDepth',...
                    'Adist',1/1.4,...
                    'Bdist',1/1,...
                    'Cdist',1/0.6,...
                    'Adeg',1.5,...
                    'Cdeg',1,...
                    'nchecks',[128 64]);

% Shrink the size of the letters so we can drop the FOV
for ii = 1:length(thisSE.recipe.assets)
    if (strcmp(thisSE.recipe.assets(ii).name,'A') || ...
       strcmp(thisSE.recipe.assets(ii).name,'B') || ...
       strcmp(thisSE.recipe.assets(ii).name,'C'))
        thisSE.recipe.assets(ii).scale = [0.5;0.5;0.5];
    end
end
%}

%%
% A small FOV is required to see the difference between the models.
thisSE.set('fov',40); 

% Set general parameters
thisSE.set('spatial resolution',[128,128]);
thisSE.set('rays per pixel',128);
thisSE.set('chromatic aberration',8);

%% Try the Navarro eye model

% This tell isetbio which model to use.
thisSE.set('model name','Navarro');

% The Navarro model has accommodation, but let's set it to infinity for now
% since other models may not have accommodation modeling.
thisSE.set('accommodation',0);

% Render!
thisSE.set('name','oi navarro'); % The name of the optical image
oiNavarro = thisSE.render();

% Show the retinal image
% Everything is very out of focus since the accommodation is set to
% infinity. 
oiWindow(oiNavarro);

% You can see the lens file used and the dispersion curves of the material
% in the working directory. However, the lens file is not written out until
% time of rendering.
% thisScene.workingDir

%% Try the Gullstrand-LeGrand Model

% The gullstrand has no accommodation modeling. 
thisSE.modelName = 'LeGrand';
thisSE.name = 'LeGrand'; % The name of the optical image
oiLeGrand = thisSE.render();
oiWindow(oiLeGrand);

%% Try Arizona eye model

thisSE.modelName = 'Arizona';
thisSE.accommodation = 0;

% Render!
thisSE.name = 'arizona'; % The name of the optical image
oiArizona = thisSE.render();
oiWindow(oiArizona);

%%
