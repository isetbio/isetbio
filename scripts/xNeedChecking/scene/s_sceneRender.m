%% s_sceneRender
% 
% Read in a multispectral scene and render it under different lights.
%
% D50.mat, D55.mat, D65.mat, D75.mat
% In the first part of this script, we use an arbitrary method for
% displaying the spectral data (map different wavelength bands into r, g
% and b)
%
% Another way to think of rendering is to calculate the xyz values for the
% spectral data. We do this by first calculating the radiance signal for
% the surfaces under daylight (think of this as color balancing) Then we
% calculate the xyz values and use an xyz to linear srgb final step should
% be to pass the linear rgb into a display gamma lut
% 
% See also: s_Scene2SampledScene and s_sceneCompress
%
% Copyright ImagEval Consultants, LLC, 2012

%% Read in the scene
wList = [400:10:700];
fullFileName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
scene = sceneFromFile(fullFileName ,'multispectral',[],[],wList);

% Have a look at the image (just mapping different spectral bands into rgb)
vcAddAndSelectObject(scene); sceneWindow;

% Plot the illuminant
plotScene(scene,'illuminant photons roi')


%% Transform the current illuminant to daylight
% notice that daylight is defined only out to ~700 nm
% try to find a spectral power distribution for daylight out to 950 nm

% Read illuminant energy.
wave  = sceneGet(scene,'wave');
daylight = ieReadSpectra('D75.mat',wave);

% Adjust function.  In this case daylight is a vector of illuminant
% energies at each wavelength.
scene = sceneAdjustIlluminant(scene,daylight);
scene = sceneSet(scene,'illuminantComment','Daylight (D75) illuminant');

% Have a look
vcAddAndSelectObject(scene); sceneWindow;
plotScene(scene,'illuminant photons roi')


%% End
