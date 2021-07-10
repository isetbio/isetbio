%% Ashby chapter
%
% cMosaic - the new cone mosaic that starts the parallel pools
% coneMosaic - The rectangular mosaic that I wrote originally
%

%% Open up the app if you want some of those images

ISETBioCSFGenerator

%%
scene = sceneFromFile;
sceneWindow(scene);

nBands = 10;
rgb = cell(nBands,1);
for ii=1:nBands
    newWave = (450:10:470) + (ii-1)*20;
    tmp = sceneSet(scene,'wave',newWave);
   rgb{ii} = sceneGet(tmp,'rgb');
end

for ii=1:nBands
    ieNewGraphWin; imshow(rgb{ii});
end

%%   Noise distribution graphs

lambda = [1 5 20];
nSamp = 10000;

% Poisson
ieNewGraphWin;

val = zeros(nSamp,numel(lambda));
for ii=1:numel(lambda)
    val(:,ii) = iePoisson(lambda(ii),'nSamp',[nSamp,1]);
end

h = histogram(val(:,1),'Normalization','pdf'); hold on;
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor = [0.8 0.8 0.8];

h = histogram(val(:,2),'Normalization','pdf'); 
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor = [0.5 0.5 0.5];

h = histogram(val(:,3),'Normalization','pdf');
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor =  [0.2 0.2 0.2];

grid on
title('Poisson');
xlabel('Counts');
ylabel('Probability density')

% Normal

ieNewGraphWin;

val = zeros(nSamp,numel(lambda));
for ii=1:numel(lambda)
    val(:,ii) = sqrt(lambda(ii))*randn(nSamp,1) + lambda(ii);
end

h = histogram(val(:,1),'Normalization','pdf'); hold on;
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor = [0.8 0.8 0.8];
h = histogram(val(:,2),'Normalization','pdf'); 
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor = [0.5 0.5 0.5];
h = histogram(val(:,3),'Normalization','pdf');
h.BinLimits = [-5 40]; h.BinWidth = 1; % h.FaceColor = [0.2 0.2 0.2];
grid on
title('Normal');
xlabel('Counts');
ylabel('Probability density')

%% Calculate cone fundamentals

% Draft of an app

theCones = coneMosaic;
theCones.wave = 400:1:700;
theLens = Lens;
theLens.wave = 400:1:700;
wave = theCones.wave;

theCones.macular.density = 0.35; % 0.35 is Default
theCones.macular.eccDensity(0:1:15)
theLens.density = 1; % 1 is default

%% 
% ieNewGraphWin;
% plot(wave,theCones.pigment.quantaFundamentals);

%% Macular pigment and photopigment
% ieNewGraphWin;
% plot(wave,theCones.qe);

%%
ieNewGraphWin([],'wide');
subplot(1,2,1)
theCones.macular.density = theCones.macular.eccDensity(0);
coneFundamentals = bsxfun(@times, theLens.transmittance', theCones.qe);
plot(wave,coneFundamentals,'Linewidth',3);
xlabel('Wavelength (nm)');
ylabel('Absorptance')
grid on;
legend('L-cones','M-cones','S-cones');
title('Central fovea');

subplot(1,2,2)
theCones.macular.density = theCones.macular.eccDensity(10); % 0.35 is Default
coneFundamentals = bsxfun(@times, theLens.transmittance', theCones.qe);
plot(wave,coneFundamentals,'Linewidth',3);
xlabel('Wavelength (nm)');
ylabel('Absorptance')
grid on;
title('10 deg eccentricity');

%% PSFs ...???

% I stopped the ISETBioCSFGenerator inside the cMosaic.oiEnsembleGenerate
% function and saved the key variables in the 'projects' directory.
%
% save('oiData','obj','oiSamplingGridDegs','varargin');

% Apparently, we can load the objects and parameters in.
chdir(fullfile(isetbioRootPath,'projects'));
foo = load('oiData');

[X,Y] = foo.obj.oiEnsembleGenerate(foo.oiSamplingGridDegs,foo.varargin{:});

% Then this runs
load('oiData','obj','oiSamplingGridDegs','varargin');
[oiEnsemble, psfEnsemble] = obj.oiEnsembleGenerate(oiSamplingGridDegs,varargin{:});

% To see the code and learn more, 
%    edit cMosaic.oiEnsembleGenerate
%
% The wavelengths and images are in the psfEnsemble variable that is
% returned.

%{
varargin =

  1×14 cell array

  Columns 1 through 6

    {'zernikeDataBase'}    {'Polans2015'}    {'subjectID'}    {[10]}    {'pupilDiameterMM'}    {[3]}

  Columns 7 through 13

    {'subtractCentral…'}    {[0]}    {'zeroCenterPSF'}    {[1]}    {'flipPSFUpsideDown'}    {[1]}    {'wavefrontSpatia…'}

  Column 14

    {[301]}
%}