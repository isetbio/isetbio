%% s_sceneFromMultispectral
%
% sceneFromFile creates a scene from multispectral and hyperspectral scene data files
%
% Copyright ImagEval Consulting LLC, 2013

%% Select the data file
%
% ISET includes a few multispectral datafiles
%
% Other data can be downloaded from 
%
% * http://imageval.com/scene-database/ 
% * http://scien.stanford.edu/index.php/hyperspectral-image-data/

fullFileName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs');

%% Read in the multispectral/hyperspectral image data
% 
wList = [400:10:700]; % You can specify the range and step size of the wavelengths you sample

scene = sceneFromFile(fullFileName,'multispectral',[],[],wList);

%% Display the scene in the scene window
vcAddAndSelectObject(scene); sceneWindow

%% Notes about multispecteal scene data files
%
% Multispectral and hyperspectral scene data files can be downloaded from
%
%  * http://imageval.com/scene-database/
%  * http://scien.stanford.edu/index.php/resources/image-databases-2/
%
% The scene data are sometimes compressed using the singular value decomposition.
%
% When the data are compressed, the scene data files contain the following information:
%
% * basis.basis: a [N x M] matrix where N is the number of wavelength samples and 
%                M is the number of spectral basis functions
% * basis.wave: a [1 x N] vector containing the wavelength samples
% * coefficients: a [rows x cols x M] matrix containing the coefficients per image pixel
% * wave: a [1 x N] vector containing the sampled wavelengths
% * illuminant.data: a [1 x N] vector describing the spectral power distribution of the illuminant 
% * illuminant.wavelength: a [1 x N] vector containing the wavelength samples
%
% When the data are uncompressed, the scene data files contain the following information: 
%
% * photons: a [rows x cols x N] matrix containing radiance (expressed in units of photons)
% * wave: a [1 x N] vector containing the sampled wavelengths
% * illuminant: a [1 x N] vector describing the spectral power distribution of the illuminant 


