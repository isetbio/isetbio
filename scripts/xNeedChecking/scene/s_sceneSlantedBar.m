% s_slantedBarScenes
%
% Create a slanted bar scene and set different scene properties. The
% slanted bar is used in measuring the ISO 12233 standard for spatial
% resolution.  This function and related ones are embedded into the
% interface and illustrated in other scripts.
%
%   sceneCreate('slantedBar',imageSize,edgeSlope, fieldOfView, illPhotons);
%
% See also:  sceneCreate, sceneSPDScale, s_adjustIlluminant,
%            s_spatialResolution
%
% Copyright ImagEval Consultants, LLC, 2010.

%% Settable slanted bar parameters
sz          = 256;
barSlope    = 2.6;
fieldOfView = 2;
meanL = 100;

%% Create an example slanted bar

% The default slanted bar is created with an illuminant of equal photons
% across wavelengths
scene = sceneCreate('slantedBar', sz, barSlope, fieldOfView);
scene = sceneAdjustLuminance(scene,meanL);

% Have a look at the image in the scene Window
vcAddAndSelectObject(scene);
sceneWindow;

% Here is the scene energy
plotScene(scene,'illuminant energy roi')

%% Change the slanted bar to a D65 illuminant, rather than equal energy

scene = sceneAdjustIlluminant(scene,'D65.mat');

% Have a look
vcAddAndSelectObject(scene);
sceneWindow;
plotScene(scene,'illuminant energy roi')

%% Create slanted bar with another slope

barSlope    = 3.6;
sz          = 128;
fieldOfView = 0.5

scene = sceneCreate('slantedBar', sz, barSlope, fieldOfView);
vcAddAndSelectObject(scene);
sceneWindow;

