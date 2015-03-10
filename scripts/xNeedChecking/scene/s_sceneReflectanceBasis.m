%% s_sceneReflectanceBasis

% Goal is to do the analysis and manipulation of the hyperspectral scene data
%   in reduced basis space

% 1.	Read in the hyperspectral scene data
% 2.	Convert the scene data to reflectance
% 3.	Find a set of basis functions for reflectance
%       Students can analyze basis functions to identify pigments, skin
%       components, and implement image enhancement algorithms in basis space
% 4.	Convert scene reflectance basis to scene radiance basis by multiplying by illuminant
% 5.	Convert radiance to xyz
% 6.	Convert xyz to rgb


% 1.	Read in the hyperspectral scene data
% fullFileName = fullfile('C:\Users\joyce\Documents\Matlab\SVN\iset-4.0\data\scenes\JoyceCloseUp.mat');
fullFileName = 'eye.mat'; % fullfile('C:\Users\joyce\Documents\Matlab\SVN\iset-4.0\data\scenes\eye.mat');
load(fullFileName);
vcAddAndSelectObject(scene); sceneWindow;
plotScene(scene,'illuminantPhotons')

%% 2.	Convert the scene data to reflectance

%% 3.	Find a set of basis functions for reflectance

%% 4.	Convert scene reflectance basis to scene radiance basis by multiplying by illuminant
% be sure to store the illuminant 

%% 5.	Convert radiance to xyz
row = sceneGet(scene,'rows');
col = sceneGet(scene,'cols');
e = Quanta2Energy(wave,double(scene.data.photons));
xyz = ieXYZFromEnergy(e,wave);
%xyz = XW2RGBFormat(xyz,row,col);

%% 6.	Convert xyz to rgb
% We find the max Y and normalize xyz when we call the function.  This is
% expected in the standard (see Wikipedia page)
Y = xyz(:,:,2); maxY = max(Y(:))
sRGB = xyz2srgb(xyz/maxY);
% Visualize the result
vcNewGraphWin; image(sRGB);
axis image;
