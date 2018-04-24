function [scene, sSamples, reflectance, rcSize] = sceneReflectanceChart(...
    sFiles, sSamples, pSize, wave, grayFlag, sampling)
% Create a reflectance chart for testing
%
% Syntax:
%	[scene, sampleList, reflectances, rcSize] = sceneReflectanceChart(...
%       sFiles, sSamples, pSize, [wave], [grayFlag=1], [sampling])
%
% Description:
%    Create a reflectance chart for testing.
%
%    If a specific set of samples is chosen they are written out in row
%    first ordering, [ 1 4 7; 2 5 8; 3 6 9]
%
%    Examples are contained in the code. To access, type 'edit
%    sceneReflectanceChart.m' into the Command Window.
%
% Inputs:
%	 sFiles       - Cell array of file names with reflectance spectra. The
%                   surfaces are drawn from the cell array named sFiles{}.
%                   It is also possible to set sFiles to a matrix comprised
%                   of reflectances.
%	 sSamples     - This can either be:
%         - A vector indicating how many surfaces to sample from each file
%         - A cell array of specifying the list of samples from each file 
%	 pSize        - (Optional) The number of pixels on the side of each
%                   square patch. Default 32.
%    wave         - (Optional) Wavelength Samples. Default pull from scene.
%    grayFlag     - (Optional) Fill the last part of the chart with gray
%                   surfaces (20% reflectance). Default 1.
%    sampling     - (Optional) Whether or not to replace samples. Default
%                   'r' (Yes to replacing by default).
%
% Outputs:
%    scene        - The reflectance chart as a scene
%    sSamples     - A cell array of the surfaces from each file, as above.
%    reflectances - The actual reflectances.
%    rcSize       - The row & column sizes.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    macbethChartCreate
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010.
%    02/01/18  jnm  Formatting
%    04/07/18  dhb  Got example to run by removing reference to a file of 
%                   skin reflectances that wasn't being found.

% Examples:
%{
    sFiles = cell(1, 3);
    sFiles{1} = fullfile(isetbioDataPath, 'surfaces', ...
        'reflectances', 'MunsellSamples_Vhrel.mat');
    sFiles{2} = fullfile(isetbioDataPath, 'surfaces', ...
        'reflectances', 'Food_Vhrel.mat');
    sFiles{3} = fullfile(isetbioDataPath, 'surfaces', ...
        'reflectances', 'DupontPaintChip_Vhrel.mat');
    sSamples = [12, 12, 25]*5; nSamples = sum(sSamples);
    pSize = 24; 

    [scene, samples] = sceneReflectanceChart(sFiles, sSamples, pSize);
    scene = sceneAdjustLuminance(scene, 100);
    vcAddAndSelectObject(scene);
    sceneWindow;
%}

if notDefined('sFiles')
    error('Surface files required'); 
else
    nFiles = length(sFiles);
end

if notDefined('pSize'), pSize = 32; end
if notDefined('grayFlag'), grayFlag = 1; end
if notDefined('sampling'), sampling = 'r'; end %With replacement by default

% Default scene
scene = sceneCreate;
if notDefined('wave')
    wave = sceneGet(scene, 'wave');
else
    scene = sceneSet(scene, 'wave', wave);
end
nWave = length(wave);
defaultLuminance = 100;  % cd/m2

% Get the reflectance samples
if iscellstr(sFiles)
    if notDefined('sSamples'), error('Surface samples required'); end
    if length(sSamples) ~= nFiles
        error('Mis-match between number of files and sample numbers');
    end
    [reflectance, sSamples] = ...
        ieReflectanceSamples(sFiles, sSamples, wave, sampling);
    
    % sSamples might be a vector, indicating the number of samples, or a
    % cell array specifying which samples.
    if iscell(sSamples)
        nSamples = 0;
        for ii = 1:nFiles, nSamples = length(sSamples{ii}) + nSamples; end
    else
        nSamples = sum(sSamples);
    end
else
    reflectance = sFiles;
    nSamples = size(reflectance, 2);
end

% Spatial arrangement
r = ceil(sqrt(nSamples));
c = ceil(nSamples /r);

% reflectance is in wave x surface format. We fill up the end of the matrix
% with gray surface reflectances.
if grayFlag
    % Create a column of gray surfaces, 20 percent reflectance
    g = 0.2 * ones(nWave, r);
    reflectance = [reflectance, g];
    c = c + 1;
end
rcSize = [r, c];

% Convert the scene reflectances into photons, assuming you have an equal
% energy illuminant.
ee = ones(nWave, 1);  % Equal energy vector
e2pFactors = Energy2Quanta(wave, ee);  % Energy to photon factor

% Illuminant
illuminantPhotons = diag(e2pFactors)*ones(nWave, 1);

% Convert the reflectances into photons
% Data from first file are in the left columns, second file next set of
% cols, and so forth. There may be a gray strip at the end.
% Scale reflectances by incorporating energy to photon scale factr
radiance = diag(e2pFactors) * reflectance;    

% Put these into the scene data structure. These are in photon units, but
% they are not scaled to reasonable photon values.
sData = reshape([radiance 0.2 * illuminantPhotons(:, ...
            ones(r * c - size(radiance, 2), 1))], [nWave r c]);
sData = permute(sData, [2 3 1]);

% Build up the size of the image regions - still reflectances
sData = imageIncreaseImageRGBSize(sData, pSize);

% Add data to scene, using equal energy illuminant
scene = sceneSet(scene, 'photons', sData);
scene = sceneSet(scene, 'illuminantPhotons', illuminantPhotons);
scene = sceneSet(scene, 'illuminantComment', 'Equal energy');
scene = sceneSet(scene, 'name', 'Reflectance Chart (EE)');
% vcAddAndSelectObject(scene); sceneWindow;

% Adjust the illuminance to a default level in cd/m2
scene = sceneAdjustLuminance(scene, defaultLuminance);
% vcAddAndSelectObject(scene); sceneWindow;

end
