% s_sceneReflectanceCharts
%
%  Create a reflectance chart with many surfaces. Such a chart can be
%  useful for testing color algorithms.  
%
% Copyright ImagEval Consultants, LLC, 2010.

%% Basic example

% Randomly select reflectances

% The files containing the reflectances are in ISET format, readable by 
% s = ieReadSpectra(sFiles{1});
sFiles = cell(1,4);
sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
sFiles{3} = fullfile(isetRootPath,'data','surfaces','reflectances','DupontPaintChip_Vhrel.mat');
sFiles{4} = fullfile(isetRootPath,'data','surfaces','reflectances','HyspexSkinReflectance.mat');

% The number of samples from each of the data sets, respectively
sSamples = [12,12,24,24];    % 

% How many row/col spatial samples in each patch (they are square)
pSize = 24;    % Patch size
wave =[];      % Whatever is in the file
grayFlag = 0;  % No gray strip
sampling = 'no replacement';
scene = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayFlag,sampling);

% Show it on the screen
vcAddAndSelectObject(scene); sceneWindow;

%% Change the illumination
% Change from the default illuminant (equal energy) to D65

wave = sceneGet(scene,'wave');  d65 = ieReadSpectra('D65',wave);

sceneD65 = sceneAdjustIlluminant(scene,d65);

sceneD65 = sceneSet(sceneD65,'name','Reflectance Chart D65');

vcAddAndSelectObject(sceneD65); sceneWindow;

%% Add a gray strip column
grayStrip = 1;

sceneGray = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayStrip);
sceneGray = sceneSet(sceneGray,'name','Reflectance Chart EE Gray Strip');

vcAddAndSelectObject(sceneGray); sceneWindow;

%% Store the parameters needed to make exactly the same chart

[sceneOriginal, storedSamples] = sceneReflectanceChart(sFiles,sSamples,pSize);
sceneOriginal = sceneSet(sceneOriginal,'name','Original');
vcAddAndSelectObject(sceneOriginal); sceneWindow;

sceneReplica = sceneReflectanceChart(sFiles,storedSamples,pSize);
sceneReplica = sceneSet(sceneReplica,'name','Replica');
vcAddAndSelectObject(sceneReplica); sceneWindow;

%% End

