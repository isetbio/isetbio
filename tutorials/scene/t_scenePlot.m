% Demonstrate the scenePlot function calls
%
% Description:
%    This tutorial is an introduction to how scenePlot may be used to
%    examine various aspects of the scene structure.
%
% See Also:
%    scenePlot
%

% History:
%    xx/xx/13  BW   Imageval Consulting, 2013
%    09/07/18  jnm  Formatting

%% Initialize ISETBio and the data
ieInit;

scene = sceneCreate;
vcAddAndSelectObject(scene);
sceneWindow;

%% Show luminance as a mesh on linear or log scale
scenePlot(scene, 'luminance mesh linear');
scenePlot(scene, 'luminance mesh log');

%% Plot spectral radiance along middle row
% Plot has horizontal position and wavelength
middleRow = round(sceneGet(scene, 'rows') / 2);
scenePlot(scene, 'hline radiance', [1, middleRow]);

%% Fourier transform of the luminance in the row
uData = scenePlot(scene, 'luminance fft hline', [1, middleRow]);

%% Radiance image with an overlaid spatial grid
gridSpacing = 21;
scenePlot(scene, 'radiance image with grid', [], gridSpacing);

%% Scene illuminant, twice
scenePlot(scene, 'illuminant photons roi');
scenePlot(scene, 'illuminant photons');

%% Depth map
uData = scenePlot(scene, 'depth map');

%% Reflectance data from an ROI
roiRect = [26 40 13 16];
uData = scenePlot(scene, 'reflectance roi', roiRect);

%% Chromaticities from within ROI
% This is a spatially simple scene, so there are just a few
scenePlot(scene, 'chromaticity', roiRect);

%% Repeat with a different ROI
roiRect = [6 51 8 12];
uData = scenePlot(scene, 'reflectance', roiRect);
scenePlot(scene, 'chromaticity', roiRect);
