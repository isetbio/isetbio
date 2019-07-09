% Create & save a set of hexagonal cone mosaics
%
% Description:
%    These have eccentricity-based cone spacing.  We save them out and put
%    them in the database for easy retrieval for demonstrations
%
%    These have an eccentricity dependent cone spacing and include an
%    S-cone free region, and a desired S-cone spacing.
%
%    Separate scripts use these mosaics to dynamic compute isomerizations
%    for different stimuli.
%
% See Also:
%   t_conesMosaixHexReg
%   advancedTutorials/t_conesMosaicHex1
%   advancedTutorials/t_conesMosaicHex6
%
%

%% Initialize
% If you want to clear all by default, or not, use
%
%   ieSessionSet('init clear',false);
%   ieSessionGet('init clear')
%
% See doc ieSessionSet() for other parameters you can control and
% ieSessionGet() for general ISETBio parameters you can retrieve.
%
ieInit;

%% Local storage location

chdir(fullfile(isetbioRootPath,'local','mosaics'));

% Flywheel set up 
st = scitran('stanfordlabs');
project = st.lookup('wandell/ISETBio Mosaics');
fovea = project.sessions.findOne('label=fovea');

%% Set mosaic FOV list.
fovList = 1;


%% Set mosaic parameters
% The various mosaic parameters and their descriptions
%
%    'name'                   - String. The name of the mosaic.
%    'resamplingFactor'       - Numeric. Sets underlying pixel spacing;
%                               controls the rectangular sampling of the
%                               hex mosaic grid.
%    'eccBasedConeDensity'    - Boolean. Whether to have an eccentricity
%                               based, spatially - varying density.
%    'sConeMinDistanceFactor  - Numeric. Min distance between neighboring
%                               S-cones = f * local cone separation - used
%                               to make the S-cone lattice semi-regular.
%    'sConeFreeRadiusMicrons' - Numeric. Radius of S-cone free retina, in
%                               microns (here set to 0.15 deg).
%    'spatialDensity'         - Vector. The KLMS vector with a LMS density
%                               of of 6:3:1.
mParams = struct(...
    'name', 'rect mosaic', ...
    'center',[3,3], ...
    'eccentricityunits','deg',...
    'spatialDensity', [0 .6 .3 .1]);

%% Create the rectangular mosaic

rectMosaic = coneMosaic(mParams);

% Set its field of view
mParams.fovdegs = 5;
rectMosaic = rectMosaic.setSizeToFOV(mParams.fovdegs);
% rectMosaic.window('show','cone mosaic');
support = [];
spread = 8;
conePlot(rectMosaic.coneLocs * 1e6, rectMosaic.pattern, support, spread);
axis off; axis image; 

%% Save the mosaic for later analysis 
mosaicFileName = sprintf('rectMosaic-center(%2.1f,%2.1f)-fov(%2.2f)',...
    mParams.center(1),mParams.center(2),mParams.fovdegs);
mosaicDataFile = sprintf('%s.mat',mosaicFileName);
mosaicMetaDataFile = sprintf('%s.json',mosaicFileName);
save(mosaicDataFile, 'rectMosaic');
jsonwrite(mosaicMetaDataFile,mParams)
[p,n,~] = fileparts(mosaicDataFile);
mosaicConesPDF = fullfile(p,[n,'-cones.pdf']);
NicePlot.exportFigToPDF(mosaicConesPDF, gcf, 300);

%% Upload to FLywheel
acqLabel = sprintf('%2.2f deg',mParams.fovdegs);
try
    acq = fovea.acquisitions.findOne(['label=',acqLabel]);
catch
    acq = fovea.addAcquisition('label',acqLabel);
end
acq.uploadFile(mosaicDataFile);
acq.uploadFile(mosaicMetaDataFile);
acq.uploadFile(mosaicConesPDF);

%% END