% s_sceneReflectanceChartBasisFunctions
%
% This script creates a color chart from a set of reflectances. The color
% chart is a scene. We find the spectral basis functions that describe
% 99.9% of the variance in the scene reflectances.
%
% Copyright Imageval Consulting LLC, 2012

%% 
s_initISET

%% Randomly select reflectances

% The files containing the reflectances are in ISET format, readable by
% s = ieReadSpectra(sFiles{1});
sFiles = cell(1,6);
sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
sFiles{3} = fullfile(isetRootPath,'data','surfaces','reflectances','DupontPaintChip_Vhrel.mat');
sFiles{4} = fullfile(isetRootPath,'data','surfaces','reflectances','HyspexSkinReflectance.mat');
sFiles{5} = fullfile(isetRootPath,'data','surfaces','reflectances','Nature_Vhrel.mat');
sFiles{6} = fullfile(isetRootPath,'data','surfaces','reflectances','Objects_Vhrel.mat');

% The number of samples from each of the data sets, respectively
sSamples = [12,12,24,5,24,12];    %

% How many row/col spatial samples in each patch (they are square)
pSize = 24;    % Patch size
wave =[];      % Whatever is in the file
grayFlag = 0;  % No gray strip
sampling = 'no replacement';
scene = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayFlag,sampling);

% Show it on the screen
vcAddAndSelectObject(scene); sceneWindow;

%% Approximate the reflectance chart with a linear model

wave        = sceneGet(scene,'wave');
reflectance = sceneGet(scene,'reflectance');

% Could subsample spatial reflectance samples if you are worried about
% memory
% reflectance = reflectance(1:3:end,1:3:end,:);

% Do not remove the mean, and require explaining 0.999 of the variance
mType = 'canonical';
bType = 0.999;
[~, basisData,~,varExplained] = hcBasis(reflectance,bType,mType);
fprintf('Variance explained %.03f by %d bases\n',...
    varExplained,size(basisData,2));
 
%% Show the basis functions
vcNewGraphWin;
plot(wave, basisData);

%% Set a lower requirement for variance explained
bType = 0.95;
[~, basisData,~,varExplained] = hcBasis( reflectance,bType,mType);
fprintf('Variance explained %.03f by %d bases\n',...
    varExplained,size(basisData,2));
vcNewGraphWin;
plot(wave, basisData);

%% Set a  requirement by number of bases
bType = 5;
[~, basisData,~,varExplained] = hcBasis( reflectance,bType,mType);
fprintf('Variance explained %.03f by %d bases\n',...
    varExplained,size(basisData,2));
tmp = size(basisData);
vcNewGraphWin;
plot(wave, basisData);

%% Could write compute and write out using sceneToFile
%
% Or, we can write out directly using
%
%  ieSaveMultiSpectralImage(fname,coef,basis,comment,imgMean,illuminant);
%

%% End

