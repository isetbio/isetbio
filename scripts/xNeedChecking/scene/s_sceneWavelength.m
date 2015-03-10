%% s_sceneWavelength
%
% Illustrate how to adjust the wavelength representation in a scene.
%
% (c) Imageval Consulting, LCC 2012

%%
s_initISET

%% 
scene = sceneCreate;
sceneGet(scene,'wave')
vcAddAndSelectObject(scene); sceneWindow;
fprintf('Note the wavelength representation in the window\n');

%% Adjust the wavelength to 5 nm spacing
scene = sceneSet(scene,'wave',400:5:700);
scene = sceneSet(scene,'name','5 nm spacing');
vcAddAndSelectObject(scene); sceneWindow;
fprintf('Note the wavelength representation in the window\n');

%%  Now get a narrow band representation
scene = sceneSet(scene,'wave',500:2:600);
scene = sceneSet(scene,'name','2 nm narrow band spacing');
vcAddAndSelectObject(scene); sceneWindow;
fprintf('Note the wavelength representation in the window\n');

%% End