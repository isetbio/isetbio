function macbethChartObject = ...
    macbethChartCreate(patchSize, patchList, spectrum)
% Initiate a scene structure of the Gretag/Macbeth Color Chart reflectances
%
% Syntax:
%	macbethChartObject = ...
%       macbethChartCreate([patchSize], [patchList], [spectrum]);
%
% Description:
%    Initiate a scene of spectral reflectance functions of the Macbeth
%    chart. The chart is coded as if you are looking at it with four rows
%    and six columns. The white surface is on the left, and the black
%    surface is on the right.
% 
%    The numbering of the surfaces starts at the upper left and counts down
%    the first column. The white patch, therefore, is number 4. The
%    achromatic (gray) series is 4:4:24.
%   
%    The surface reflectance function information are contained in the slot
%    macbethChartObject.data. The data file contains spectral information
%    between 380 to 1068 nm.
%
%    There are examples contained in the code. To access the examples, type
%    'edit macbethChartCreate.m' into the Command Window.
%
% Inputs:
%    patchSize          - (Optional) The patch size of the mcc color
%                         sections. Default is 16 pixels.
%    patchList          - (Optional) The list of the color patches. Default
%                         is 1:24 (The list of all macbeth surfaces).
%    spectrum           - (Optional) The color spectrum. The default is to
%                         initialize a default (hyperspectral) spectrum.
%                         This can be the wavelengths of the colors
%                         (400:10:700).
%
% Outputs:
%	 macbethChartObject - The created mcc object
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/01/18  jnm  Formatting

% Examples:
%{
    patchSize = 16;
    patchList = 1:24;
    macbethChartObject = macbethChartCreate(patchSize, patchList);

    % To read a different spectral range
    spectrum.wave = (370:10:730);
    macbethChartObject = macbethChartCreate([], [], spectrum);

    % To make an image
    macbethChartObject = macbethChartCreate(patchSize, patchList);
    spd = macbethChartObject.data;
    imageSPD(spd);
%}
%{
    % To read the chart reflectances 
    patchSize = 1;
    patchList = 1:24;
    macbethChartObject = macbethChartCreate(patchSize, patchList);
    r = macbethChartObject.data;
    reflectances = RGB2XWFormat(r)';
    vcNewGraphWin;
    plot(reflectances)
%}
%{
    % To read just the gray series
    patchSize = 16;
    patchList = 4:4:24;
    macbethChartObject = macbethChartCreate(patchSize, patchList);
    spd = macbethChartObject.data;
    imageSPD(spd);
%}

%% Initialize object
% The object is basically a scene, though it is missing some fields.
% These get filled in by sceneCreate, which uses this routine.
macbethChartObject.name = 'Macbeth Chart';
macbethChartObject.type = 'scene';

% This is the size in pixels of each Macbeth patch
if notDefined('patchSize'), patchSize = 16; end

% These are the patches we are trying to get
% If we want just the gray series we can set patchList = 19:24;
if notDefined('patchList'), patchList = 1:24; end

%% Surface reflectance spectrum
if notDefined('spectrum')
    macbethChartObject = initDefaultSpectrum(...
        macbethChartObject, 'hyperspectral');
else
    macbethChartObject = sceneSet(macbethChartObject, ...
        'spectrum', spectrum);
end

% Read wavelength information from the macbeth chart data
wave = sceneGet(macbethChartObject, 'wave');
nWaves = sceneGet(macbethChartObject, 'nwave');

% Read the MCC reflectance data
fName = fullfile(isetbioDataPath, 'surfaces', 'macbethChart.mat');
macbethChart = ieReadSpectra(fName, wave);

% Sort out whether we have the right set of patches
macbethChart = macbethChart(:, patchList);
if length(patchList) == 24
    macbethChart = reshape(transpose(macbethChart), 4, 6, nWaves);
else
    macbethChart = reshape(transpose(macbethChart), 1, ...
        length(patchList), nWaves);
end

% Make it the right patch size
macbethChartObject.data = imageIncreaseImageRGBSize(...
    macbethChart, patchSize);

end