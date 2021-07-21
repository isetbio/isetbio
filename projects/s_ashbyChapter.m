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

%% Cone fundamentals

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

%% PSFs using the ISETBIOCsfGenerator methods

% I stopped the ISETBioCSFGenerator inside the 
%
%     cMosaic.oiEnsembleGenerate
%
% method.  While there, I saved the key optics variables in the 'projects'
% directory. 
%
% save('oiData','obj','oiSamplingGridDegs','varargin');

% We can load the objects and parameters.
chdir(fullfile(isetbioRootPath,'projects'));

% Load the variables
load('oiData','obj','oiSamplingGridDegs','varargin');

% The key variable is obj, the ISETBio cMosaic
%
% obj.oiEnsembleGenerate runs to create the psfs given the parameters
%
% edit cMosaic.oiEnsembleGenerate 
%
%  brings you to the code
%
% The meaning of the parameters is:
%    centerPos = oiSamplingGridDegs;
%{
varargin =
    {'zernikeDataBase'}    {'Polans2015'}    
    {'subjectID'}            {[10]}    
    {'pupilDiameterMM'}      {[3]}    
    {'subtractCentral…'}     {[0]}    
    {'zeroCenterPSF'}        {[1]}    
    {'flipPSFUpsideDown'}    {[1]}    
    {'wavefrontSpatia…'}    {[301]}
% Convert the varargin to a struct that we can edit
args = struct;
for ii=1:2:numel(varargin)
    args.(varargin{ii}) = varargin{ii+1};
end
%}

% Choose Artalright eye, subject 1?? 74, 26.  or 85
centerPos = [0,1];

args.zernikeDataBase = 'Polans2015';

args.zernikeDataBase = 'Artal2012';
args.subjectID       = 26;
args.pupilDiameterMM = 3;
args.subtractCentralRefraction = false;
args.zeroCenterPSF = true;
args.flipPSFUpsideDown = true;
% args.deNoisedZernikeCoefficients = false;
args.wavefrontSpatialSamples = 301;

[oiEnsemble, psfEnsemble] = obj.oiEnsembleGenerate(centerPos,args);
thisPSF = psfEnsemble{1};

ieNewGraphWin;  
thisW = 7;   % Select a wavelength, % 3 is 450 nm, 7 is 550 nm 
h = mesh(thisPSF.supportX,thisPSF.supportY,thisPSF.data(:,:,thisW)/max2(thisPSF.data(:,:,thisW)));
tMarks  = [-15:5:15]; mn = -5; mx = 5;
set(gca,'xtick',tMarks,'xlim',[mn mx],'ylim',[mn mx]);
set(gca,'ytick',tMarks);
xlabel('Pos (microns)');
ylabel('Pos (microns)');
zlabel('Relative intensity');

title(sprintf('Wave %d Sub %d',thisPSF.supportWavelength(thisW),args.subjectID));

%% To get cone diameter as a function of eccentricity

ecc = (0.1:0.2:20);
sz  = 


%% Make an additive Gaussian noise and stimulus dependent noise

% scene = sceneCreate('linear intensity ramp');
scene = sceneCreate('lstar');

scene = sceneSet(scene,'fov',5);

oi = oiCreate;
oi = oiCompute(oi,scene);

%%  Photon noise
eTime = 0.002;
cDensity = [0 0 1 0];
fov = 6;

cmP = coneMosaic;
cmP.setSizeToFOV(fov);
% All M cones
cmP.spatialDensity = cDensity;
cmP.integrationTime = eTime;   % 10 ms
cmP.compute(oi);

photonNoise = cmP.absorptions;
mx = max(photonNoise(:));
ieNewGraphWin; imagesc(photonNoise,[0 mx]); 
axis image; colormap(hot(128));
title('Poisson');brighten(0.3);

cmP.name = 'photon';
% cmP.window;

%%
cmG = coneMosaic;
cmG.setSizeToFOV(fov);
cmG.spatialDensity = cDensity;
cmG.integrationTime = eTime;  
cmG.noiseFlag = 'none';
cmG.compute(oi);

noNoise = cmG.absorptions;

mn = max(noNoise(:));
equivNoise = randn(size(noNoise))*sqrt(mn);

gaussNoise = noNoise + equivNoise;
gaussNoise(gaussNoise<0) = 0;
mx = max(gaussNoise(:));

ieNewGraphWin; imagesc(gaussNoise,[0 mx]); 
axis image; colormap(hot(128));
title('Gaussian'); brighten(0.3);

cmG.absorptions = gaussNoise;
cmG.name = 'gaussian';
% cmG.window;
