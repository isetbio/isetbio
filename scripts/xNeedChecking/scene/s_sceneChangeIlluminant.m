%% s_sceneChangeIlluminant
%
% Illustrate how to  adjust the illuminant of the current scene, simulating
% a change in the spectral power distribution.
%
% You can set the illuminant to one of the standard SPDs in the directory
% data/lights using the GUI pulldown
%
%   Edit | Adjust SPD | Change illuminant
%
% See also: s_Exercise, sceneAdjustIlluminant,s_illuminantCorrection
%
% Copyright ImagEval Consultants, LLC, 2010.

%%
s_initISET

%% Create a scene
%  
% sceneCreate is a function that creates a scene
% If there are no arguments, sceneCreate will simulate a Macbeth
% ColorChecker uniformly illuminated with daylight (D65)

scene = sceneCreate;

% Have a look at the image
vcAddAndSelectObject(scene); sceneWindow;

% Plot the illuminant
plotScene(scene,'illuminant photons roi')


%% Replace the current scene illuminant with Tungsten

% Read illuminant energy.
wave  = sceneGet(scene,'wave');
TungstenEnergy = ieReadSpectra('Tungsten.mat',wave);

% Adjust function.  In this case TungstenEnergy is a vector of illuminant
% energies at each wavelength.
scene = sceneAdjustIlluminant(scene,TungstenEnergy);
scene = sceneSet(scene,'illuminantComment','Tungsten illuminant');

% Have a look
vcAddAndSelectObject(scene); sceneWindow;
plotScene(scene,'illuminant photons roi')

%% Create a scene based on multispectral scene data
%
%  see s_sceneFromMultispectral.m
sceneFile = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
scene = sceneFromFile(sceneFile,'multispectral');
scene = sceneAdjustLuminance(scene,61); % This sets the mean scene luminance
scene = sceneSet(scene,'fov',26.5); % match the scene field of view (fov) with the sensor fov

vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow
plotScene(scene,'illuminant energy roi')

%% Set illuminant to equal energy

scene = sceneAdjustIlluminant(scene,'equalEnergy.mat');

vcAddAndSelectObject('scene',scene); sceneWindow; % display sceneWindow

%% Convert the scene to the sunset color, Horizon_Gretag
%
% Notice that in this case 'Horizon_Gretag.mat' is a file name, not a
% data vector. 

scene = sceneAdjustIlluminant(scene,'Horizon_Gretag.mat');
vcAddAndSelectObject('scene',scene); sceneWindow; 


%% End
