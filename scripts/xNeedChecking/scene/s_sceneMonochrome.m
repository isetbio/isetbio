%% s_sceneMonochrome
%
% Read in a monochrome image and convert it to a unispectral image.
% Specify the display as a crt
%
% 
% Copyright ImagEval Consultants, LLC, 2005.

dispFile = 'crt';
d = displayCreate(dispFile);

fName = 'cameraman.tif';
scene = sceneFromFile(fName,'monochrome',100,dispFile);

% Adjust the unispectral to D65
bb = blackbody(sceneGet(scene,'wave'),6500,'energy');
scene = sceneAdjustIlluminant(scene,bb);

vcAddAndSelectObject(scene)
sceneWindow


%%