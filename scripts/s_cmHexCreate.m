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

fovList = [2,4];

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
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...
    'eccBasedConeDensity', true, ...
    'eccBasedConeQuantalEfficiency', true, ...
    'eccBasedMacularPigment', false, ...
    'sConeMinDistanceFactor', 3.0, ...
    'sConeFreeRadiusMicrons', 0.15 * 300, ...
    'spatialDensity', [0 .6 .3 .1]);

% larger than default tolerances to speed-up computation. For production
% work, either do not set, or set to equal or lower than 0.01.
mParams.quality.tolerance1 = 0.5;
% larger than default tolerances to speed-up computation, For production
% work, either do not set, or set to equal or lower than 0.001.
mParams.quality.tolerance2 = 0.05;
% How much larger lattice to generate so as to minimize artifacts in cone
% spacing near the edges. If empty, a dynamic adjustment of margin is done
% for mosaics < 1.0 degs.
mParams.quality.marginF = [];
% Iterations
mParams.quality.gridAdjustmentIterations = 50;
% Some day, this will be a meaningful parameter
mParams.center = 'center(0,0)';

for pIndex = 1:numel(fovList)
    mParams.fovDegs = fovList(pIndex);
    mosaicFileName = sprintf('hexMosaic-%s-fov(%2.2f)', mParams.center,mParams.fovDegs);
    
    tic
    %% Generate the mosaic. This takes a little while.
    hexMosaic = coneMosaicHex(mParams.resamplingFactor, ...
        'name', mParams.name, ...
        'fovDegs', mParams.fovDegs, ...
        'eccBasedConeDensity', mParams.eccBasedConeDensity, ...
        'eccBasedConeQuantalEfficiency', mParams.eccBasedConeQuantalEfficiency, ...
        'eccBasedMacularPigment', mParams.eccBasedMacularPigment, ...
        'sConeMinDistanceFactor', mParams.sConeMinDistanceFactor, ...
        'sConeFreeRadiusMicrons', mParams.sConeFreeRadiusMicrons, ...
        'spatialDensity', mParams.spatialDensity, ...
        'latticeAdjustmentPositionalToleranceF', mParams.quality.tolerance1, ...
        'latticeAdjustmentDelaunayToleranceF', mParams.quality.tolerance2, ...
        'maxGridAdjustmentIterations', mParams.quality.gridAdjustmentIterations, ...
        'marginF', mParams.quality.marginF);
    toc
    
    % Save the mosaic for later analysis (Why are we saving with this
    % version?)
    mosaicDataFile = sprintf('%s.mat',mosaicFileName);
    mosaicMetaDataFile = sprintf('%s.json',mosaicFileName);
    save(mosaicDataFile, 'hexMosaic');
    jsonwrite(mosaicMetaDataFile,mParams)
    
    
    %% Print mosaic info
    % hexMosaic.displayInfo();
    
    %% Visualize the mosaic, showing inner segment and geometric area.
    % The inner segment being the light collecting area
    
    % Choose aperture from 'both', 'lightCollectingArea', 'geometricArea'
    visualizedAperture = 'lightCollectingArea';
    hFig1 = hexMosaic.visualizeGrid(...
        'axesHandle', gca, ...
        'visualizedConeAperture', visualizedAperture, ...
        'apertureShape', 'disks', ...
        'ticksInMicrons', true);
    set(hFig1, 'Position', [50 50 1200 1200]);

    [p,n,~] = fileparts(mosaicDataFile);
    mosaicConesPDF = fullfile(p,[n,'-cones.pdf']);
    NicePlot.exportFigToPDF(mosaicConesPDF, hFig1, 300);

    %% Visualize mosaic w/ overlayed theoretical & measured cone dens plots
    
    % coneDensityContour levels are in cones/mm^2

    % maxEccMicrons = 300*mParams.fovDegs;
    contourLevels = 1000 * linspace(20,200,5);
    hFig2 = hexMosaic.visualizeGrid(...
        'axesHandle', gca, ...
        'visualizedConeAperture', visualizedAperture, ...
        'apertureShape', 'disks', ...
        'labelConeTypes', false, ...
        'overlayHexMesh', true, ...
        'overlayConeDensityContour', 'theoretical_and_measured', ...
        'coneDensityContourLevels', contourLevels, ...
        'ticksInMicrons', true);
    set(hFig2, 'Position', [50 50 1200 1200]);

    [p,n,e] = fileparts(mosaicDataFile);
    mosaicDensityPDF = fullfile(p,[n,'-density.pdf']);
    NicePlot.exportFigToPDF(mosaicDensityPDF, hFig2, 300);

    %% Upload to FLywheel
    acqLabel = sprintf('%2.2f deg',mParams.fovDegs);
    try
        acq = fovea.acquisitions.findOne(['label=',acqLabel]);
    catch
        acq = fovea.addAcquisition('label',acqLabel);
    end
    acq.uploadFile(mosaicDataFile);
    acq.uploadFile(mosaicMetaDataFile);
    acq.uploadFile(mosaicConesPDF);
    acq.uploadFile(mosaicDensityPDF);
    
end % pIndex

%% END